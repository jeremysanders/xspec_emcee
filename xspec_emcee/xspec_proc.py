from __future__ import print_function, division, absolute_import

import atexit
import os.path
import os
import re
import subprocess

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
            cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            universal_newlines=True, bufsize=1,
        )

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
        self.buffer += os.read(self.popen.stdout.fileno(), 8192).decode('utf8')
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
