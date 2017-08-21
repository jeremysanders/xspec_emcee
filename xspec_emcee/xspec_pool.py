from __future__ import print_function, division, absolute_import

import itertools
import select
from collections import defaultdict

import numpy as N

class CombinedModel:
    """Model containing all xspec models to evaluate to give a
    total model."""

    def __init__(self, xspecmodels):
        self.xspecmodels = xspecmodels
        self.update_thawed()

    def update_thawed(self):
        existing = set()
        self.thawedparams = []
        for model in self.xspecmodels:
            for tp in model.thawedparams:
                if tp not in existing:
                    self.thawedparams.append(tp)
                    existing.add(tp)

    def log_norms_priors(self, minnorm=1e-10):
        """Modify priors on norms to be flat in log space."""

        def getPrior(par):
            def prior(v):
                if v < par.minval or v > par.maxval:
                    return -N.inf
                return -N.log(v)
            return prior

        for par in self.thawedparams:
            if par.name == 'norm':
                print(' Using prior for log value on parameter %i:%s:%s:%i' % (
                        par.xspecindex, par.model, par.cmpt, par.index))

                par.minval = max(par.minval, minnorm)
                par.prior = getPrior(par)

    def prior(self, vals):
        """Calculate total prior for parameters."""
        return sum((par.prior(val)
                    for par, val in itertools.izip(self.thawedparams, vals)))

    def updateParams(self, vals):
        """Update thawed parameters with model parameters."""
        for par, val in itertools.izip(self.thawedparams, vals):
            par.currentval = val

    def linkParameters(self, linkexpr):
        """Link two parameters.

        Form is
        [xcmindex1:][[modelname1]:]paramindex1 = [xcmindex2:][[modelname2]:]paramindex2

        Default xcmindex is 1 and default modelname is unnamed
        """

        def defpart(t):
            p = t.split(':')
            if len(p) > 3:
                raise RuntimeError('Parameter expression should have at most 3 parts')
            elif len(p) == 2:
                p = [1] + p
            elif len(p) == 1:
                p = [1, 'unnamed'] + p
            p[0] = int(p[0])
            p[1] = p[1].strip() if p[1].strip() else 'unnamed'
            p[2] = int(p[2])
            return p

        def findindex(params, modname, paramindex):
            for i, p in enumerate(params):
                if p.model == modname and p.index == paramindex:
                    return i
            raise RuntimeError('Cannot find parameter', modname, paramindex)

        left, right = linkexpr.split('=')
        left, right = defpart(left), defpart(right)

        print(" Linked parameter %i:%s:%i to %i:%s:%i" % tuple(right+left))

        # get xspec models
        lxmodel = self.xspecmodels[left[0]-1]
        rxmodel = self.xspecmodels[right[0]-1]

        # get parameter indices
        lthaw, rthaw = lxmodel.thawedparams, rxmodel.thawedparams
        lidx = findindex(lthaw, left[1], left[2])
        ridx = findindex(rthaw, right[1], right[2])
        lp, rp = lthaw[lidx], rthaw[lidx]
        print("  (%s:%s:%s:%s -> %s:%s:%s:%s)" % (
                right[0], rp.model, rp.cmpt, rp.name,
                left[0], lp.model, lp.cmpt, lp.name))

        # assign right parameter to left model
        lthaw[lidx] = rthaw[ridx]

        self.update_thawed()

class ProcState:
    """This object is for handling the processing state for a set of
    XspecModel objects."""

    def __init__(self, combmodel, xmodel, paramlist, likes, toprocess):
        self.combmodel = combmodel
        self.xmodel = xmodel
        self.paramlist = paramlist
        self.likes = likes
        self.toprocess = list(toprocess)

        # map fileno to xspec process
        self.fileno_to_proc = {x.fileno(): x for x in xmodel.procs}

        # fileno which are free to process
        self.free = list(self.fileno_to_proc.keys())

        # filenos which are doing work
        self.processing = {}

    def _check(self):
        """Check whether processes have returned results and get them."""

        # iterate over processes which have written to stdout
        for fileno in select.select(list(self.processing.keys()), [], [])[0]:
            proc = self.fileno_to_proc[fileno]
            result = proc.read_buffer()
            if result is not None:
                # valid result, so get likelihood
                like = -0.5*float(result)
                self.likes[self.processing[fileno]] += like

                # free up process for next job
                self.free.append(fileno)
                del self.processing[fileno]

    def _send_job(self):
        """Send the next job to process."""
        paridx = self.toprocess.pop()
        paramset = self.paramlist[paridx]

        self.combmodel.updateParams(paramset)

        # build up newpar command to send to xspec
        modparams = defaultdict(list)
        for param in self.xmodel.thawedparams:
            mpm = modparams[param.model]
            while len(mpm) < param.index-1:
                mpm.append('')
            mpm.append('%e' % param.currentval)
        # newpar command for each model
        cmds = []
        for model, pars in modparams.iteritems():
            cmd = 'newpar %s1-%i & %s' % (
                '' if model == 'unnamed' else model+':',
                len(pars), ' & '.join(pars))
            cmds.append(cmd)
        # command to get output statistic
        cmds.append('emcee_tcloutr stat')

        fileno = self.free.pop()
        proc = self.fileno_to_proc[fileno]
        proc.send_cmd('\n'.join(cmds))
        self.processing[fileno] = paridx

    def loop_iter(self):
        """Does a cycle of starting new jobs and getting the results
        of old ones.

        Returns None if none left."""

        if self.free and self.toprocess:
            self._send_job()

        elif self.processing:
            # check for a completed process
            self._check()

        return not self.toprocess and not self.processing

class XspecPool:
    def __init__(self, combmodel):
        """Fake pool object to return likelihoods for parameter sets."""

        self.combmodel = combmodel

        # keep track of evaluations
        self.itercount = 0

    def map(self, dummyfunc, paramlist):
        """Return a list of lnprob values for the list parameter sets
        given.

        Note: dummyfunc is never called!
        """

        # get prior for initial likelihood
        likes = N.array([self.combmodel.prior(v) for v in paramlist])
        # list of parameters with finite priors
        toprocess = list(N.nonzero(N.isfinite(likes))[0])

        #notfinite = ~N.isfinite(likes)
        #for i in N.nonzero(notfinite)[0]:
        #    print(paramlist[i])

        # each xspecmodel has a processing state for each of the parameters
        states = [
            ProcState(self.combmodel, xmodel, paramlist, likes, toprocess)
            for xmodel in self.combmodel.xspecmodels
            ]

        while True:
            alldone = True
            for state in states:
                done = state.loop_iter()
                alldone = alldone and done
            if alldone:
                break

        likefilt = likes[N.isfinite(likes)]
        if len(likefilt) > 0 and self.itercount % 2 == 0:
            print('%5i   mean=<%9.1f> max=<%9.1f> std=<%9.1f> good=<%4i/%4i>' % (
                    self.itercount // 2,
                    likefilt.mean(),
                    likefilt.max(),
                    likefilt.std(),
                    len(likefilt), len(likes),
                    ))
        self.itercount += 1

        return likes
