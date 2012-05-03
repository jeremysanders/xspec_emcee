# Module for controlling remote xspec process
# Jeremy Sanders 2012

import subprocess
import os
import select
import atexit
import re
import sys
import os.path
import socket
import select

thisdir = os.path.dirname(os.path.abspath(__file__))

def deleteFile(filename):
    """Delete file, ignoring errors."""
    try:
        os.unlink(filename)
    except OSError:
        pass

# keep track of xspec invocations which need finishing
_finishatexit = []

class XControl(object):
    """Object controls xspec instance to get fit statistic."""

    dologging = False
    specialcode = '@@EMCEE@@'

    def __init__(self, initxcm, system=None):
        """Start xspec process and load xcm file."""

        if system is None or system == 'localhost':
            # use local system
            cmd = ['xspec']
        else:
            # run on remote system
            cmd = ['ssh', '-x', system,
                   "source /etc/profile; do_xray_profile; xspec"]

        # talk to remote xspec
        self.xspecsub = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE)

        self.poll = select.poll()
        self.poll.register(self.xspecsub.stdout.fileno(), select.POLLIN)

        self.throwAwayOutput()
        _finishatexit.append(self)

        # debug logging
        if self.dologging:
            logfile = '/home/jss/xspec.log.%i.%s' % (
                os.getpid(), system)
            deleteFile(logfile)
            self.xspecsub.stdin.write('log %s\n' % logfile)

        # load helper functions
        self.xspecsub.stdin.write('source %s\n' %
                                  os.path.join(thisdir, 'emcee_helpers.tcl'))
        # load xcm file
        self.xspecsub.stdin.write('cd %s\n' %  os.path.dirname(initxcm))
        self.xspecsub.stdin.write('@%s\n' %  os.path.basename(initxcm))

        self.buffer = ''

    def throwAwayOutput(self):
        """Ignore output from program until no more data available."""
        while True:
            if self.poll.poll(0):
                t = os.read(self.xspecsub.stdout.fileno(), 1024)
            else:
                break

    def readBlock(self):
        """Read a block of text written by command.
        Returns list of lines written between two specialcode lines
        """
        started = False
        lines = []
        while True:
            line = self.xspecsub.stdout.readline().strip()
            if line == self.specialcode:
                if started:
                    break
                else:
                    started = True
            elif started:
                lines.append(line)
        return lines

    def finish(self):
        """Exit xspec and wait for exiting."""
        self.xspecsub.stdin.write('quit\n')
        self.throwAwayOutput()
        self.xspecsub.stdout.close()
        self.xspecsub.stdin.close()
        self.xspecsub.wait()
        del _finishatexit[ _finishatexit.index(self) ]

    def getParams(self):
        """Get list of variable parameters."""

        self.throwAwayOutput()
        self.xspecsub.stdin.write('emcee_interrogate_params\n')
        lines = iter( self.readBlock() )

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
        parlist = []
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

        return parlist

    def enterLoop(self):
        """Enter continuous loop until exit"""
        self.xspecsub.stdin.write('emcee_loop\n')

    def setParameters(self, paridxs, params):
        """Set parameters in xspec and wait until statistic."""
        args = ['%i %g' % x for x in zip(paridxs, params)]
        self.xspecsub.stdin.write( ' '.join(args) + '\n' )

    def pollStatistic(self):
        """See whether there is a statistic result."""
 
        if self.poll.poll(0):
            t = os.read(self.xspecsub.stdout.fileno(), 4096)
            if not t:  # file was closed
                return None
        else:
            return None

        # add to text buffer
        self.buffer += t

        # check for a match to the buffer
        m = re.search( r'%s\s*([0-9.]+)\s*%s' % (self.specialcode,
                                                 self.specialcode),
                       self.buffer )
        if m:
            self.buffer = ''
            return float(m.group(1))
        else:
            return None

def _finishXSpecs():
    """Finish any remaining xspecs if finish() does not get called above."""
    while _finishatexit:
        _finishatexit[0].finish()

atexit.register(_finishXSpecs)
