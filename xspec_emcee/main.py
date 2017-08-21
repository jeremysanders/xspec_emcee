#!/usr/bin/env python

"""
Use EMCEE to do MCMC in Xspec.
Jeremy Sanders 2012-2017

Requires Python 2.7+, numpy, h5py and emcee
"""

from __future__ import print_function, division, absolute_import

import sys
import argparse
import time
import re
import itertools

import h5py
import numpy as N
import emcee

from .xspec_model import XspecModel
from .xspec_pool import XspecPool, CombinedModel

def gen_initial_parameters(parameters, nwalkers):
    """Construct list of initial parameter values for each walker."""
    p0 = []
    for walker in xrange(nwalkers):
        pwalker = []
        # for each walker, use initial parameters based on parameter
        # and delta parameter
        for par in parameters:
            width = par.delta
            swidth = par.sigma*0.1
            if swidth > 0 and swidth < width:
                # use sigma if delta is badly adjusted
                width = swidth

            for i in xrange(10000):
                v = N.random.normal(par.initval, width)
                if N.isfinite(par.prior(v)):
                    pwalker.append(v)
                    break
            else:
                raise RuntimeError(
                    'Could not generate initial parameter with finite prior')

        p0.append(N.array(pwalker))
    return N.array(p0)

def expand_systems(systems):
    """Allow system*N syntax in systems."""
    out = []
    for s in systems:
        m = re.match(r'([A-Za-z0-9]+)\*([0-9]+)', s)
        if m:
            out += [m.group(1)]*int(m.group(2))
        else:
            out.append(s)
    return out

def do_mcmc(xcms,
            nwalkers=100, nburn=100, niters=1000,
            systems=['localhost'],
            outchain=['out.chain'],
            outhdf5='out.hdf5',
            debug=False,
            continuerun=False,
            autosave=True,
            nochdir=False,
            initialparameters=None,
            lognorm=False,
            link=[]):
    """Do the actual MCMC process."""

    print("Loading XCM file(s)")
    xmodels = []
    for i, xcm in enumerate(xcms):
        print(' Loading', xcm)
        xmodels.append( XspecModel(
                xcm, expand_systems(systems), debug=debug, nochdir=nochdir,
                xspecindex=i+1 ))
    combmodel = CombinedModel(xmodels)

    if lognorm:
        print("Using prior equivalent to log parameter")
        combmodel.log_norms_priors()

    if link:
        print("Linking parameters")
        for expr in link:
            combmodel.linkParameters(expr)

    print("Number of free parameters: %i\n" % len(combmodel.thawedparams))

    if not initialparameters:
        print("Generating initial parameters")
        p0 = gen_initial_parameters(combmodel.thawedparams, nwalkers)
    else:
        print("Loading initial parameters from", initialparameters)
        p0 = N.loadtxt(initialparameters)

    ndims = p0.shape[1]
    pool = XspecPool(combmodel)

    # sample the mcmc
    sampler = emcee.EnsembleSampler(nwalkers, ndims, None, pool=pool)

    print("Starting MCMC")
    if not continuerun and nburn > 0:
        # burn in
        print("Burn in period started")
        pos, prob, state = sampler.run_mcmc(p0, nburn)
        sampler.reset()
        print("Burn in period finished")
    else:
        # no burn in
        state = None
        pos = p0

    if not continuerun:
        # create new datasets, extensible along number of iterations
        hdf5file = h5py.File(outhdf5, "w")
        chain = hdf5file.create_dataset(
            "chain",
            (nwalkers, niters, ndims),
            maxshape=(nwalkers, None, ndims))
        lnprob = hdf5file.create_dataset(
            "lnprob",
            (nwalkers, niters),
            maxshape=(nwalkers, None))
        start = 0

    else:
        print("Continuing from existing chain in", outhdf5)

        hdf5file = h5py.File(outhdf5, "r+")
        chain = hdf5file["chain"]
        lnprob = hdf5file["lnprob"]

        start = chain.attrs["count"]
        pos = N.array(chain[:, start-1, :])
        print("Restarting at iteration", start)

        chain.resize((nwalkers, niters, ndims))
        lnprob.resize((nwalkers, niters))

    # iterator interface allows us to trap ctrl+c and know where we are
    lastsave = time.time()
    index = start
    try:
        for p, l, s in sampler.sample(
            pos, rstate0=state, storechain=False,
            iterations=niters-start):

            chain[:, index, :] = p
            lnprob[:, index] = l
            index += 1

            if autosave and time.time() - lastsave > 60*10:
                chain.attrs["count"] = index
                hdf5file.flush()

    except KeyboardInterrupt:
        chain.attrs["count"] = index
        print("Ctrl+C pressed - ending")

    else:
        chain.attrs["count"] = index
        write_xspec_chains(outchain, chain, lnprob, combmodel)

    hdf5file.close()

def write_xspec_chains(filenames, chain, lnprob, combmodel):
    """Write an xspec text chain file for each xcm input file."""

    nwalkers, niters, ndims = chain.shape
    # real length could be shorter
    niters = chain.attrs["count"]

    def innerwrite(chainf, xmodel):
        chainf.write('! Markov chain file generated by xspec "chain" command.\n')
        chainf.write('!    Do not modify, else file may not reload properly.\n')

        chainf.write('!Length: %i  Width: %i\n' % (
                niters*nwalkers, len(xmodel.thawedparams)+1))
        chainf.write('!Type: GoodmanWeare\n')
        chainf.write('!NWalkers: %i\n' % nwalkers)

        # header for contents of file
        hdr = []
        for par, idx in itertools.izip(
            xmodel.thawedparams, xmodel.xspec_thawed_idxs()):
            hdr.append("%s %s %s" % (
                    idx, par.name,
                    par.unit if par.unit else "0"))
        hdr.append("Likelihood")
        chainf.write('!%s\n' % ' '.join(hdr))

        fmt = '\t'.join(['%g']*len(xmodel.thawedparams))

        # write each walker separately
        for wi in xrange(nwalkers):
            chainw = chain[wi, :, :]
            statw = lnprob[wi, :]

            # then the iterations
            for params, stat in itertools.izip(chainw, statw):
                combmodel.updateParams(params)
                vals = [p.currentval for p in xmodel.thawedparams]
                line = fmt % tuple(vals) + '\t' + '%g' % stat + '\n'
                chainf.write(line)

    for chainfilename, xmodel in zip(filenames, combmodel.xspecmodels):
        print('Writing text chain', chainfilename)
        with open(chainfilename, 'w') as chainf:
            innerwrite(chainf, xmodel)

def run():
    """Main program."""

    p = argparse.ArgumentParser(
        description="Xspec MCMC with EMCEE. Jeremy Sanders 2012-2017.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument("xcms", metavar="XCM", nargs="+",
                   help="Input XCM file")
    p.add_argument("--niters", metavar="N", type=int, default=5000,
                   help="Number of iterations")
    p.add_argument("--nburn", metavar="N", type=int, default=500,
                   help="Number of burn iterations")
    p.add_argument("--nwalkers", metavar="N", type=int, default=50,
                   help="Number of walkers")
    p.add_argument("--systems", default="localhost", metavar="LIST",
                   help="Space-separated list of computers to run on")
    p.add_argument("--output-hdf5", default="emcee.hdf5", metavar="FILE",
                   help="Output HDF5 file")
    p.add_argument("--output-chain",
                   metavar="FILE", action="append",
                   help="Output text file")
    p.add_argument("--continue-run",  action="store_true", default=False,
                   help="Continue from an existing chain (in HDF5)")
    p.add_argument("--debug", action="store_true", default=False,
                   help="Create Xspec log files")
    p.add_argument("--no-chdir", action="store_true", default=False,
                   help="Do not chdir to XCM file directory before execution")
    p.add_argument("--initial-parameters", metavar="FILE",
                   help="Provide initial parameters")
    p.add_argument("--log-norm", action="store_true", default=False,
                   help="Use priors equivalent to using log norms")
    p.add_argument('--chunk-size', metavar='N', type=int, default=4,
                   help='Currently ignored')
    p.add_argument("--link", metavar="EXPR", action="append",
                   help="Link two parameters in model")

    args = p.parse_args()

    # get list of output chain files
    # this is complex as it may not be the same as the number of xcm files
    outchain = args.output_chain
    if outchain is None:
        if len(args.xcms) == 1:
            outchain = ['emcee.chain']
        else:
            outchain = ['emcee.chain.%i' % (i+1) for i in range(len(args.xcms))]
    else:
        if len(outchain) == 1 and '%' in outchain[0]:
            outchain = [outchain[0] % (i+1) for i in range(len(args.xcms))]
        else:
            if len(outchain) != len(args.xcms):
                raise RuntimeError('Requires same number of output chains as input chains')

    sampler = do_mcmc(
        args.xcms,
        systems = args.systems.split(),
        nwalkers = args.nwalkers,
        nburn = args.nburn,
        niters = args.niters,
        outchain = outchain,
        outhdf5 = args.output_hdf5,
        continuerun = args.continue_run,
        debug = args.debug,
        nochdir = args.no_chdir,
        initialparameters = args.initial_parameters,
        lognorm = args.log_norm,
        link = args.link,
    )

    print("Done")

if __name__ == '__main__':
    run()
