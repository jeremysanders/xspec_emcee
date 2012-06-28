import os
import os.path
import re
import itertools
import time

import numpy as N

import subprocessing

# script to start xspec
thisdir = os.path.dirname( os.path.abspath(__file__) )
start_xspec = os.path.join(thisdir, 'start_xspec.sh')
# helper routines to load
xspec_helpers = os.path.join(thisdir, 'emcee_helpers.tcl') 

# get result (stripped of whitespace)
result_re = re.compile(r'@@EMCEE@@\s*(.*?)\s*@@EMCEE@@', flags=re.DOTALL)

class XspecPool(subprocessing.Pool):
    def __init__(self, xcm, systems, debug=False):
        cmds = [ [start_xspec, system] for system in systems ]
        self.xcm = os.path.abspath(xcm)
        self.debug = debug

        # start up processes
        subprocessing.Pool.__init__(self, cmds)

        # get parameters for xcm model
        self._get_params()

        # enter execution loop
        for popen in self.popens:
            popen.stdin.write('emcee_loop\n')

        # feedback
        self.count = 0

    def _get_params(self):
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
            if len(vals) == 1:
                # switch parameter
                vals = [vals[0], -1, -1e99, -1e99, 1e99, 1e99]

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
                   'cmpt': cmpts[i+1]}
            parlist.append(par)

        # identify unlinked and unfrozen parameters
        self.paridxs = []   # indices
        self.parvals = []   # values
        self.parameters = []   # full details of parameter
        for par in parlist:
            if not par['linked'] and par['val_delta'] > 0.:
                self.paridxs.append(par['index'])
                self.parvals.append(par['val_init'])
                self.parameters.append(par)

    def init_subprocess(self, popen):
        """Initialise the xspec process given."""

        # load helper routines
        popen.stdin.write('source %s\n' % xspec_helpers)

        if self.debug:
            popen.stdin.write('set logfile '
                              '$env(HOME)/xspec.log.[info hostname].[pid]\n')
            popen.stdin.write('file delete $logfile\n')
            popen.stdin.write('log $logfile\n')
            popen.stdin.write('set EMCEE_DEBUG 1\n')

        # load xcm in current directory
        popen.stdin.write('cd %s\n' % os.path.dirname(self.xcm))
        popen.stdin.write('@%s\n' % os.path.basename(self.xcm))

        time.sleep(0.5)

    def close_subprocess(self, popen):
        """Tell loop to exit."""
        popen.stdin.write('quit\n')
        subprocessing.Pool.close_subprocess(self, popen)

    def identify_lnprob(self, text):
        """Is the result in the output?"""
        match = result_re.search(text)
        if match:
            statistic = float( match.group(1) )
            if statistic < 0:
                return -N.inf
            else:
                return -statistic * 0.5
        return None

    def in_bounds(self, params):
        """Are the parameters in bounds?"""

        # return no probability if out of bounds
        for pidx, pval in itertools.izip(self.paridxs, params):
            xp = self.parlist[pidx-1]
            if pval < xp['val_hardmin'] or pval > xp['val_hardmax']:
                return False
        return True

    def send_parameters(self, stdin, params):
        """Send parameters to process."""

        if self.count % 100 == 0:
            print self.count
        self.count += 1

        if not self.in_bounds(params):
            # tell process to retn error to ourselves!
            stdin.write('returnerror\n')
        else:
            # set the parameters and await result
            maxidx = max(*self.paridxs)
            cmd = [ "1-%i" % maxidx ]

            paridx = 0
            pairs = zip(self.paridxs, params)
            pairs.sort()
            for i in xrange(maxidx):
                if i == pairs[paridx][0]-1:
                    cmd.append(str(pairs[paridx][1]))
                    paridx += 1
                else:
                    cmd.append('')

            text = ' & '.join(cmd)
            stdin.write( text + '\n' )
