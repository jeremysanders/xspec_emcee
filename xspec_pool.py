from __future__ import print_function, division

import os
import os.path
import re
import itertools
import time
import atexit
import subprocess
import select

import numpy as N

# script to start xspec
thisdir = os.path.dirname( os.path.abspath(__file__) )
start_xspec = os.path.join(thisdir, 'start_xspec.sh')
# helper routines to load
xspec_helpers = os.path.join(thisdir, 'emcee_helpers.tcl') 

# get result (stripped of whitespace)
result_re = re.compile(r'@@EMCEE@@\s*(.*?)\s*@@EMCEE@@', flags=re.DOTALL)

class XspecPool(object):
    def __init__(self, xcm, systems, debug=False, nochdir=False,
                 lognorm=False, chunksize=4):
        cmds = [ [start_xspec, system] for system in systems ]
        self.xcm = os.path.abspath(xcm)
        self.debug = debug
        self.nochdir = nochdir
        self.lognorm = lognorm
        self.chunksize = chunksize

        # list of open subprocesses
        self.popens = []
        self.buffers = {}
        for cmd in [ [start_xspec, system] for system in systems ]:
            self.popens.append(self.init_subprocess(cmd))
            self.buffers[self.popens[-1]] = ''
        atexit.register(self.finish)

        # get parameters for xcm model
        self.get_params()

        # enter execution loop
        for popen in self.popens:
            popen.stdin.write('emcee_loop\n')

        # feedback
        self.count = 0

    def finish(self):
        """Finish all processes."""
        # tell processes to finish
        for p in self.popens:
            p.stdin.write('quit\n')
            p.stdin.close()
        # wait until they have closed
        for p in self.popens:
            p.wait()
        del self.popens[:]

    def get_params(self):
        """Get list of parameters from first process."""
        p = self.popens[0]

        # send command to get parameters
        p.stdin.write('emcee_interrogate_params\n')
        txt = ''
        while True:
            txt += os.read(p.stdout.fileno(), 4096)
            match = result_re.search(txt)
            if match:
                break
        lines = iter( match.group(1).split('\n') )

        # build up mapping from parameter index -> component name
        cmpts = {}
        ncmpt = int(lines.next())
        for i in xrange(ncmpt):
            cmpinfo = lines.next().split()
            name = "%s[%i]" % (cmpinfo[0], i+1)
            start, count = int(cmpinfo[1]), int(cmpinfo[2])
            for c in xrange(start, start+count):
                cmpts[c] = name

        # now build up a list of parameter dicts
        self.parlist = parlist = []
        nvars = int(lines.next())
        for i in xrange(nvars):
            name = lines.next().split()
            link = lines.next().split()
            vals = [float(x) for x in lines.next().split()]
            sigma = float(lines.next())
            log = name==['norm'] and self.lognorm

            if log:
                # convert to log
                for p in 0, 2, 3, 4, 5:
                    v = vals[p]
                    if v<=0:
                        v = 1e-99
                    vals[p] = N.log10(v)

            par = {'index': i+1,
                   'name': name[0],
                   'unit': '' if len(name) == 1 else name[1],
                   'linked': link[0].lower() == 't',
                   'val_init': vals[0],
                   'val_delta': vals[1],
                   'val_hardmin': vals[2],
                   'val_softmin': vals[3],
                   'val_softmax': vals[4],
                   'val_hardmax': vals[5],
                   'val_sigma': sigma,
                   'cmpt': cmpts[i+1],
                   'log': log,
                   }
            parlist.append(par)

        # identify unlinked and unfrozen parameters
        self.paridxs = []   # indices
        self.parvals = []   # values
        self.parameters = []   # full details of parameter
        self.hardmin = []
        self.hardmax = []
        self.parlog = []
        for par in parlist:
            if not par['linked'] and par['val_delta'] > 0.:
                self.paridxs.append(par['index'])
                self.parvals.append(par['val_init'])
                self.hardmin.append(par['val_hardmin'])
                self.hardmax.append(par['val_hardmax'])
                self.parlog.append(par['log'])
                self.parameters.append(par)
        self.hardmin = N.array(self.hardmin)
        self.hardmax = N.array(self.hardmax)
        self.parlog = N.array(self.parlog, dtype=N.bool)
        self.maxidx = max(*self.paridxs)

    def init_subprocess(self, cmd):
        """Initialise the xspec process given."""

        popen = subprocess.Popen(
            cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

        # load helper routines
        popen.stdin.write('source %s\n' % xspec_helpers)

        if self.debug:
            popen.stdin.write('set logfile '
                              '$env(HOME)/xspec.log.[info hostname].[pid]\n')
            popen.stdin.write('file delete $logfile\n')
            popen.stdin.write('log $logfile\n')
            popen.stdin.write('set EMCEE_DEBUG 1\n')

        # load xcm in current directory
        if self.nochdir:
            popen.stdin.write('@%s\n' % self.xcm)
        else:
            popen.stdin.write('cd %s\n' % os.path.dirname(self.xcm))
            popen.stdin.write('@%s\n' % os.path.basename(self.xcm))

        time.sleep(0.5)
        return popen

    def get_parameters(self, params):
        """Send parameters to process."""
        if N.any(params < self.hardmin) or N.any(params > self.hardmax):
            return None
        else:
            # set the parameters and await result
            cmd = ["1-%i" % self.maxidx] + [""]*self.maxidx
            paramscpy = N.array(params)
            paramscpy[self.parlog] = 10**paramscpy[self.parlog]

            for i, p in itertools.izip(self.paridxs, paramscpy):
                cmd[i] = str(p)
            return " & ".join(cmd)

    def map(self, function, paramlist):
        """Return a list of lnprob values for the list parameter sets
        given.

        Note: function is never called!
        """

        # convert parameters to text
        retn = N.zeros(len(paramlist), dtype=N.float64)*N.nan
        numbers = []
        params = []
        for i, p in enumerate(paramlist):
            ptext = self.get_parameters(p)
            if ptext is None:
                retn[i] = -1
            else:
                params.append(ptext)
                numbers.append(i)

        chunksize = self.chunksize

        free = range(len(self.popens))
        waiting = {}

        while len(params) > 0 or len(waiting) > 0:
            # send off parameters to any waiting processes
            while free and len(params)>0:
                idx = free.pop()
                popen = self.popens[idx]
                chunk = params[:chunksize]
                chunkparams = ' '.join('{%s}' % p for p in chunk)
                popen.stdin.write('batch %s\n' % chunkparams)
                waiting[popen.stdout.fileno()] = (idx, popen, numbers[:chunksize])
                params = params[chunksize:]
                numbers = numbers[chunksize:]

            # have any finished processing?
            stdouts = select.select(waiting.keys(), [], [])[0]
            for stdout in stdouts:
                idx, popen, pnumbers = waiting[stdout]
                self.buffers[popen] += os.read(stdout, 4096)
                match = result_re.search(self.buffers[popen])
                if match:
                    results = match.group(1)
                    results = N.fromstring(results, dtype=N.float64, sep=' ')
                    retn[pnumbers] = results
                    self.buffers[popen] = ''
                    free.append(idx)
                    del waiting[stdout]

        assert N.all(N.isfinite(retn))

        resfilt = retn[retn > 0]
        if len(resfilt) > 0 and self.count % 2 == 0:
            print('%5i   mean=<%9.1f> min=<%9.1f> bad=<%4i/%4i>' % (
                    self.count // 2, resfilt.mean(), resfilt.min(),
                    N.sum((retn<0).astype(N.intc)), len(retn),
                    ))
        self.count += 1

        lnprob = retn*-0.5
        lnprob[lnprob>0] = -N.inf

        return lnprob
