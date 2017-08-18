from __future__ import print_function, division

import re
import numpy as N

from xspec_proc import XspecProc

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

class XspecModel:
    """Handle multiple Xspec processes and model."""

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
