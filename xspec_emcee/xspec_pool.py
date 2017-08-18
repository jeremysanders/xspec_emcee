from __future__ import print_function, division, absolute_import

import itertools
import select
from collections import defaultdict

import numpy as N

class XspecPool(object):
    def __init__(self, xspecmodel):
        """Initialise pool.

        xspecmodel is a XspecModel object.
        """

        self.xmodel = xspecmodel

        # keep track of evaluations
        self.itercount = 0

    def map(self, dummyfunc, paramlist):
        """Return a list of lnprob values for the list parameter sets
        given.

        Note: dummyfunc is never called!
        """

        fileno_to_proc = {x.fileno(): x for x in self.xmodel.procs}

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
                for par, val in itertools.izip(self.xmodel.thawedpars, parset):
                    prior += par.prior(val)
                if not N.isfinite(prior):
                    likeout[paridx] = -N.inf
                    continue
                likeout[paridx] = prior

                # build up newpar command to send to xspec
                modpars = defaultdict(list)
                for par, val in itertools.izip(self.xmodel.thawedpars, parset):
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
