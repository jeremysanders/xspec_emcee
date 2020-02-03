from __future__ import print_function, division, absolute_import

import re
import numpy as N

from .xspec_proc import XspecProc

class Par:
    """Model parameter convenience class."""

    def __init__(self, **argsv):
        self.__dict__.update(argsv)

        if 'prior' not in argsv:
            self.prior = self._flatPrior

    def __repr__(self):
        out = []
        for k, v in sorted(self.__dict__.items()):
            if k != 'prior':
                out.append('%s=%s' % (k, repr(v)))
        return '<Par(%s)>' % ', '.join(out)

    def _flatPrior(self, val):
        """Calculate prior log likelihood for parameter."""
        if val < self.minval or val > self.maxval:
            return -N.inf
        return 0.

class XspecModel:
    """Handle multiple Xspec processes and model."""

    def __init__(self, xcm, systems, debug=False, nochdir=False, xspecindex=-1, nofit=False):

        self.nofit = nofit
        self.xspecindex = xspecindex
        self.procs = [
            XspecProc(xcm, system, debug=debug, nochdir=nochdir)
            for system in systems
            ]
        self.models, self.pars = self._get_pars()

        # filter thawed parameters
        self.thawedparams = []
        # thawed parameter indices in xspec format
        for modelname in self.models:
            thawed = [p for p in self.pars[modelname] if p.thawed]
            self.thawedparams += thawed

    def xspec_thawed_idxs(self):
        """Return list of thawed parameter indices in xspec format."""
        return [
            ('' if p.model=='unnamed' else p.model+':')+str(p.index)
            for p in self.thawedparams
            ]

    def finish(self):
        """Finish all processes."""
        for proc in self.procs:
            self.send_finish()
        for proc in self.procs:
            self.wait_finish()
        del self.procs[:]

    def _get_pars(self):
        """Get parameters from xcm file

        Returns a list of models and a dict mapping model names to a
        list of Par objects. 'unnamed' is the main xspec model.
        """

        self.procs[0].wait()
        if not self.nofit:
            print('Initial fit...')
            # initial fit to get sigma values
            self.procs[0].send_cmd('fit')
            self.procs[0].wait()

        models = []
        modelpars = {}

        print('Interrogating model...')
        pars = self.procs[0].single_cmd('emcee_pars')
        modname = None
        for line in pars.split('\n'):
            # look for line "Model name:...Active/On"
            m = re.match(r'^\s*Model ([0-9A-Za-z_]+):.+Active/On$', line)
            if m:
                modname = m.group(1)
                models.append(modname)
                modelpars[modname] = []
                continue

            # look for line "Model ...Active/On" (default unnamed model)
            m = re.match(r'^\s*Model .+Active/On$', line)
            if m:
                modname = 'unnamed'
                models.append(modname)
                modelpars[modname] = []
                continue

            # look for parameter 1 2 name
            m = re.match(r'^\s*([0-9]+)\s+([0-9]+)\s+([A-Za-z0-9_]+).*$', line)
            if m:
                par = self._handle_par(modname, int(m.group(1)), m.group(3))
                modelpars[modname].append(par)
                continue

        if not models or not modelpars:
            raise RuntimeError('Could not find model in xcm file')

        print("Done.")
        print("Obtained %i model(s):" % len(models))
        for mod in models:
            numpars = len(modelpars[mod])
            numthawed = len([p for p in modelpars[mod] if p.thawed])
            print("  Model '%s', %i parameter(s), %i thawed" % (mod, numpars, numthawed))

        return models, modelpars

    def _handle_par(self, modname, paridx, cmptname):
        """Take modelname, parameter index and component name, and return parameter Par.,"""

        p0 = self.procs[0]
        mprefix = '' if modname=='unnamed' else modname+':'
        parinfo = p0.tclout('pinfo %s%i' % (mprefix, paridx)).split()
        linked = p0.tclout('plink %s%i' % (mprefix, paridx))[:1] == 'T'

        # parameter value, range, etc
        pvals = [
            float(x) for x in
            p0.tclout('param %s%i' % (mprefix, paridx)).split() ]

        if len(pvals) == 1:
            # switch parameter
            minval = maxval = delta = None
            thawed = False
        else:
            thawed = pvals[1] > 0 and not linked
            minval, maxval, delta = pvals[2], pvals[5], pvals[1]

        if thawed:
            sigma = float(p0.tclout('sigma %s%i' % (mprefix, paridx)))
        else:
            sigma = 0.

        par = Par(
            name=parinfo[0],
            unit='' if len(parinfo)==1 else parinfo[1],
            cmpt=cmptname,
            model=modname,
            index=paridx,
            initval=pvals[0],
            minval=minval,
            maxval=maxval,
            linked=linked,
            thawed=thawed,
            delta=delta,
            sigma=sigma,
            currentval=None,
            xspecindex=self.xspecindex,
        )
        return par

