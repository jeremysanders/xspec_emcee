# a pool class which talks to external xspecs

import time
import itertools

import numpy as N
from emcee_xcontrol import XControl

class Pool(object):

    def __init__(self, initxcm, systems):
        """Initialise the pool of remote system."""

        # initialse functions
        self.controls = []
        self.freecontrols = set()
        self.waitingcontrols = {}
        for system in systems:
            c = XControl(initxcm, system=system)
            self.controls.append(c)
            self.freecontrols.add(c)

        # get parameters and bounds
        self.params = self.controls[0].getParams()

        # identify unlinked and unfrozen parameters
        self.paridxs = []   # indices
        self.parvals = []   # values
        for par in self.params:
            if not par['linked'] and par['val_delta'] > 0.:
                self.paridxs.append(par['index'])
                self.parvals.append(par['val_init'])

        for cntrl in self.controls:
            cntrl.enterLoop()

        self.evalcount = 0

    def incCount(self):
        self.evalcount += 1
        if self.evalcount % 20 == 0:
            print self.evalcount

    def inBounds(self, params):
        """Are the parameters in bounds?"""

        # return no probability if out of bounds
        for pidx, pval in itertools.izip(self.paridxs, params):
            xp = self.params[pidx-1]
            if pval < xp['val_hardmin'] or pval > xp['val_hardmax']:
                return False
        return True

    def map(self, func, paramlist):
        """This is called by emcee to process the set of inputs.

        This code ignores the value of func!
        """

        # keep track of returned statistics
        retnresults = [None]*len(paramlist)

        # copy of index to input and parameter result
        parcopy = []
        for i, p in enumerate(paramlist):
            parcopy.append( (i, p) )

        while len(parcopy) > 0 or len(self.waitingcontrols) > 0:

            # start new jobs if possible
            while len(self.freecontrols) > 0 and len(parcopy) > 0:

                num, params = parcopy[0]
                if self.inBounds(params):
                    control = iter(self.freecontrols).next()

                    # send control new parameters
                    control.setParameters(self.paridxs, params)

                    # keep track of which result we're waiting for
                    self.waitingcontrols[control] = num
                    # no longer free
                    self.freecontrols.remove(control)
                else:
                    retnresults[num] = -N.inf
                    self.incCount()

                del parcopy[0]

            # poll running controls for returned results
            for control, num in list(self.waitingcontrols.items()):
                statistic = control.pollStatistic()
                if statistic is not None:
                    # save the result (probability is -chi^2/2)
                    retnresults[num] = -statistic/2.
                    self.incCount()

                    # we're now free for the next job
                    del self.waitingcontrols[control]
                    self.freecontrols.add(control)

            time.sleep(0.001)

        assert None not in retnresults
        return retnresults
