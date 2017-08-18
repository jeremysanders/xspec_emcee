from __future__ import print_function, division

import os
import os.path
import re
import itertools
import time
import atexit
import subprocess
import select
from collections import defaultdict

import numpy as N

# script to start xspec
thisdir = os.path.dirname( os.path.abspath(__file__) )
start_xspec = os.path.join(thisdir, 'start_xspec.sh')
# helper routines to load in xspec
xspec_helpers = os.path.join(thisdir, 'emcee_helpers.tcl') 

# get result (stripped of whitespace)
result_re = re.compile(r'>EMCEE>\s*(.*?)\s*<EMCEE<', flags=re.DOTALL)

# keep track of running xspec processes to make sure they are ended
running_procs = set()
@atexit.register
def _finish_running_procs():
    """End any xspec processes in event of crash."""
    for p in running_procs:
        p.send_finish()
    for p in list(running_procs):
        p.wait_finish()

class Par:
    """Model parameter convenience class."""

    def __init__(self, **argsv):
        self.__dict__.update(argsv)

        if 'prior' not in argsv:
            self.prior = self._flatPrior

    def __repr__(self):
        out = []
        for k, v in sorted(self.__dict__.items()):
            out.append('%s=%s' % (k, repr(v)))
        return '<Par(%s)>' % ', '.join(out)

    def _flatPrior(self, val):
        """Calculate prior log likelihood for parameter."""
        if val < self.minval or val > self.maxval:
            return -N.inf
        return 0.

class XspecProc:
    """Handle Xspec process."""

    def __init__(self, xcm, system, debug=False, nochdir=False):
        self.popen = self._init_subprocess(xcm, system, debug, nochdir)
        self.buffer = ''
        running_procs.add(self)

    def fileno(self):
        """Get fileno to wait for output."""
        return self.popen.stdout.fileno()

    def _init_subprocess(self, xcm, system, debug, nochdir):
        """Initialise the xspec process given."""

        cmd = [start_xspec, system]
        popen = subprocess.Popen(
            cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

        # load helper routines
        popen.stdin.write('source %s\n' % xspec_helpers)

        if debug:
            popen.stdin.write(
                'set logfile '
                '$env(HOME)/xspec.log.[info hostname].[pid]\n')
            popen.stdin.write('file delete $logfile\n')
            popen.stdin.write('log $logfile\n')
            popen.stdin.write('set EMCEE_DEBUG 1\n')

        # load xcm in current directory
        absxcm = os.path.abspath(xcm)
        if nochdir:
            popen.stdin.write('@%s\n' % absxcm)
        else:
            popen.stdin.write('cd %s\n' % os.path.dirname(absxcm))
            popen.stdin.write('@%s\n' % os.path.basename(absxcm))

        popen.stdin.write('emcee_startup\n')
        if not debug:
            popen.stdin.write('emcee_loop\n')
        return popen

    def send_cmd(self, cmd):
        """Send a command."""
        self.popen.stdin.write(cmd + '\n')
        self.popen.stdin.flush()

    def read_buffer(self):
        """Read from process into buffer.

        If there is a result in the buffer, then return string value
        """
        self.buffer += os.read(self.popen.stdout.fileno(), 8192)
        match = result_re.search(self.buffer)
        if match:
            result = match.group(1)
            self.buffer = ''  # assume only one result in read!
            return result
        else:
            return None

    def single_cmd(self, cmd):
        """Send command and return result."""
        self.send_cmd(cmd)
        while True:
            result = self.read_buffer()
            if result is not None:
                break
        return result

    def tclout(self, args):
        """Shortcut to get tclout results."""
        return self.single_cmd('emcee_tcloutr ' + args)

    def send_finish(self):
        """Tell subprocess to finish."""
        self.send_cmd('quit')
        self.popen.stdin.close()

    def wait_finish(self):
        """Wait for subprocess to finish."""
        self.popen.wait()
        self.popen = None
        running_procs.remove(self)

class Xspec:
    """Handle the multiple processes in xspec."""

    def __init__(self, xcm, systems, debug=False, nochdir=False):

        self.procs = [
            XspecProc(xcm, system, debug=debug, nochdir=nochdir)
            for system in systems
            ]
        self.models, self.pars = self._get_pars()

        # filter thawed parameters
        self.thawedpars = []
        # thawed parameter indices in xspec format
        for modelname in self.models:
            thawed = [p for p in self.pars[modelname] if p.thawed]
            self.thawedpars += thawed

    def xspec_thawed_idxs(self):
        """Return list of thawed parameter indices in xspec format."""
        return [
            ('' if p.model=='unnamed' else p.model+':')+str(p.index)
            for p in self.thawedpars
            ]

    def finish(self):
        """Finish all processes."""
        for proc in self.procs:
            self.send_finish()
        for proc in self.procs:
            self.wait_finish()
        del self.procs[:]

    def log_norms_priors(self, minnorm=1e-10):
        """Modify priors on norms to be flat in log space."""

        def getPrior(par):
            def prior(v):
                if v < par.minval or v > par.maxval:
                    return -N.inf
                return -N.log(v)
            return prior

        for par in self.thawedpars:
            if par.name == 'norm':
                print(' Using prior for log value on parameter %s:%s:%i' % (
                        par.model, par.cmpt, par.index))

                par.minval = max(par.minval, minnorm)
                par.prior = getPrior(par)

    def _get_pars(self):
        """Get parameters from xcm file

        Returns a list of models and a dict mapping model names to a
        list of Par objects. 'unnamed' is the main xspec model.
        """

        # initial fit required to get sigma values
        self.procs[0].send_cmd('fit')

        xmodel = self.procs[0].tclout('model')
        models = ['unnamed'] + re.findall('([A-Za-z0-9]+):', xmodel)

        # get model component information
        model_pars = {}
        for model in models:
            model_pars[model] = self._get_model_pars(model)

        return models, model_pars

    def _get_model_pars(self, model):
        """Get parameters for model given.

        xspec supports multiple models
        model is unnamed or model name
        """

        p0 = self.procs[0]
        bmodel = '' if model=='unnamed' else model
        ncmpt = int(p0.tclout('modcomp %s' % bmodel))
        ndgrp = int(p0.tclout('datagrp'))

        cmptmodelidx = 1
        pars = []
        # iterate over data groups
        for dgrp in range(ndgrp):
            # iterate over components in data group
            for cmpt in range(ncmpt):
                pars += self._get_cmpt_pars(model, dgrp, cmpt, cmptmodelidx)
                cmptmodelidx += 1
        return pars

    def _get_cmpt_pars(self, model, dgrp, cmpt, cmptmodelidx):
        """Get parameters associated with model, datagroup and component.

        cmptmodelidx is a numerical index which goes into the component
        name, counting inside the model."""

        p0 = self.procs[0]
        cmodel = '' if model=='unnamed' else model+':'
        cmptinfo = p0.tclout(
            'compinfo %s%i %i' % (cmodel, cmpt+1, dgrp+1)).split()
        startpar, npars = int(cmptinfo[1]), int(cmptinfo[2])
        cmptname = '%s<%i>' % (cmptinfo[0], cmptmodelidx)

        pars = []
        # iterate over parameters in component
        for paridx in range(startpar, startpar+npars):
            # name and unit
            parinfo = p0.tclout('pinfo %s%i' % (cmodel, paridx)).split()

            # is model linked?
            linked = p0.tclout('plink %s%i' % (cmodel, paridx))[:1] == 'T'

            # parameter value, range, etc
            pvals = [
                float(x) for x in
                p0.tclout('param %s%i' % (cmodel, paridx)).split() ]

            if len(pvals) == 1:
                # switch parameter
                minval = maxval = delta = None
                thawed = False
            else:
                thawed = pvals[1] > 0 and not linked
                minval, maxval, delta = pvals[2], pvals[5], pvals[1]

            if thawed:
                sigma = float(p0.tclout('sigma %s%i' % (cmodel, paridx)))
            else:
                sigma = 0.

            par = Par(
                name=parinfo[0],
                unit='' if len(parinfo)==1 else parinfo[1],
                cmpt=cmptname,
                model=model,
                index=paridx,
                initval=pvals[0],
                minval=minval,
                maxval=maxval,
                linked=linked,
                thawed=thawed,
                delta=delta,
                sigma=sigma,
                )
            pars.append(par)
        return pars

class XspecPool(object):
    def __init__(self, xspec):
        """Initialise pool.

        xspec is a Xspec object.
        """

        self.xspec = xspec

        # keep track of evaluations
        self.itercount = 0

    def map(self, dummyfunc, paramlist):
        """Return a list of lnprob values for the list parameter sets
        given.

        Note: dummyfunc is never called!
        """

        fileno_to_proc = {x.fileno(): x for x in self.xspec.procs}

        # file numbers of processes doing nothing
        free = list(fileno_to_proc.keys())
        # map file numbers of processes working on data to output index
        processing = {}

        # indices of parameters to be processed
        toprocess = range(len(paramlist))

        # output likelihoods
        likeout = N.zeros(len(paramlist), dtype=N.float64)

        def check():
            for fileno in select.select(list(processing.keys()), [], [])[0]:
                proc = fileno_to_proc[fileno]
                result = proc.read_buffer()
                if result is not None:
                    like = -0.5*float(result)
                    # add likelihood to prior
                    likeout[processing[fileno]] += like
                    # free up process for next job
                    free.append(fileno)
                    del processing[fileno]

        while toprocess:
            if free:
                paridx = toprocess.pop()
                parset = paramlist[paridx]

                # calculate prior for set of parameters
                # skip doing evaluation if prior is not finite
                prior = 0.
                for par, val in itertools.izip(self.xspec.thawedpars, parset):
                    prior += par.prior(val)
                if not N.isfinite(prior):
                    likeout[paridx] = -N.inf
                    continue
                likeout[paridx] = prior

                # build up newpar command to send to xspec
                modpars = defaultdict(list)
                for par, val in itertools.izip(self.xspec.thawedpars, parset):
                    mpm = modpars[par.model]
                    while len(mpm) < par.index-1:
                        mpm.append('')
                    mpm.append('%e' % val)
                # newpar command for each model
                cmds = []
                for model, pars in modpars.iteritems():
                    cmd = 'newpar %s1-%i & %s' % (
                        '' if model == 'unnamed' else model+':',
                        len(pars), ' & '.join(pars))
                    cmds.append(cmd)
                # command to get output statistic
                cmds.append('emcee_tcloutr stat')

                fileno = free.pop()
                proc = fileno_to_proc[fileno]
                proc.send_cmd('\n'.join(cmds))
                processing[fileno] = paridx

            else:
                # check for a completed process
                check()

        # wait for final completion
        while processing:
            check()

        resfilt = likeout[N.isfinite(likeout)]
        if len(resfilt) > 0 and self.itercount % 2 == 0:
            print('%5i   mean=<%9.1f> max=<%9.1f> std=<%9.1f> good=<%4i/%4i>' % (
                    self.itercount // 2,
                    resfilt.mean(),
                    resfilt.max(),
                    resfilt.std(),
                    len(resfilt), len(likeout),
                    ))
        self.itercount += 1

        return likeout
