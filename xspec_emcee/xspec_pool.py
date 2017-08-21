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
        self.thawedparams = []
        for model in xspecmodels:
            self.thawedparams += model.thawedparams

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

        # first get priors
        priors = N.array([self.combmodel.prior(v) for v in paramlist])

        # this mask is those parameters which have a finite prior
        mask = N.isfinite(priors)

        # add model likelihoods to priors
        likes = priors
        for xmodel in self.combmodel.xspecmodels:
            innerlikes = self._mapinner(xmodel, mask, paramlist)
            likes += innerlikes

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

    def _mapinner(self, xmodel, mask, paramlist):

        fileno_to_proc = {x.fileno(): x for x in xmodel.procs}

        # file numbers of processes doing nothing
        free = list(fileno_to_proc.keys())
        # map file numbers of processes working on data to output index
        processing = {}

        # indices of parameters to be processed
        toprocess = list(N.nonzero(mask)[0])

        # output likelihoods
        likes = N.zeros(len(paramlist), dtype=N.float64)

        def check():
            for fileno in select.select(list(processing.keys()), [], [])[0]:
                proc = fileno_to_proc[fileno]
                result = proc.read_buffer()
                if result is not None:
                    # return likelihood
                    like = -0.5*float(result)
                    likes[processing[fileno]] = like
                    # free up process for next job
                    free.append(fileno)
                    del processing[fileno]

        while toprocess:
            if free:
                paridx = toprocess.pop()
                paramset = paramlist[paridx]

                self.combmodel.updateParams(paramset)

                # build up newpar command to send to xspec
                modparams = defaultdict(list)
                for param in xmodel.thawedparams:
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

        return likes
