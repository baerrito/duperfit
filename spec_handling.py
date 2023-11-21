import os
import sys
import fileinput
import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage import median_filter
from auxiliary import *
from bspline_fitting import *
from placed import *

def estimate_eflux(wav,flux,lines=None,zguess=0.,silence=False,method="bspl"):
    """
    FUNCTION: ESTIMATE_EFLUX
    ****************************************************************************
    Estimates an error spectrum by computing the residuals from a given fit and
    interpolating the running standard deviation of those residuals across the
    entire spectrum.
    ****************************************************************************
    INPUT:
        wav     :   Wavelength in Angstroms | real
        flux    :   Spectral flux density | real
    OPTIONAL INPUT:
        lines   :   Line information | numpy structure, default None
            lines['ION']       --> Ion name | character
            lines['FIT']       --> True if line is ignored in a fit | bool
            lines['LINE']      --> Wavelength of line in Angstroms | real
            lines['EW']        --> Equivalent width of the line (-1 if
                                       absorption, +1 if emission) | real
            lines['TIEDTO']    --> Something | real
            lines['TIEDRATIO'] --> Something else | real
            lines['ISABS']     --> True if line is absorption | bool
            lines['MASKRAD']   --> Width for line masking in Angstroms | real
        zguess  :   Estimate of supernova redshift | real, default 0.0
        silence :   Set True to silence messages | bool, default False
        method  :   Uncertainty estimation method | character, default `bspl`
            `sgfilt`    --> Uses Savitzky-Golay filtering
            `bspl`      --> uses a B-Spline iterative fit (recommended)
    OUTPUT:
        eflux   :   Estimated error spectrum | real
    ****************************************************************************
    Written by Michael Baer
        -`bspl` method adapted from a routine by Robert Quimby
        -`sgfilt` method inspired by `bspl` method
    """
    if method == "sgfilt":
        norm=1.00/np.median(flux)
        spc=[]
        for ll,ff in zip(wav,flux):
            spc.append(np.array((ll,ff*norm),dtype=[('wav',float),('flux',float)]))
        spc=np.array(spc)
        smspc=sgfilter(spc,wav,order=3)
        # End hacks; set any NaN values to nearest
        lstore=[];lflag=0;ustore=[]
        for i,fl in enumerate(smspc['flux']):
            if np.isnan(fl) and (lflag==0):
                lstore.append(i)
            if np.isfinite(fl) and (lflag==0):
                lflag=1
            if np.isnan(fl) and (lflag==1):
                ustore.append(i)
        lstore=np.array(lstore);ustore=np.array(ustore)
        smspc['flux'][lstore]=smspc['flux'][lstore.max()+1]
        smspc['flux'][ustore]=smspc['flux'][ustore.min()-1]
        smooth=smspc['flux']
        res=flux-smooth/norm
        dwav=np.median(wav[1:]-wav[:-1])
        step=dwav/30.
        testwavs=[];stds=[]
        for thiswav in np.arange(np.min(wav),np.max(wav)-step*0.5,step):
            w=(wav>=thiswav)&(wav<thiswav+step)
            nw=np.sum(w)
            if nw<10:
                continue
            sres=res[w]
            sres=np.sort(sres)
            std=np.std(sres[int(0.1*(nw-1)):int(0.9*(nw-1))])
            medres=np.median(sres)
            w2=np.where(np.abs(sres-medres)<3.0*std)
            testwavs.append(np.mean(wav[w]))
            stds.append(np.std(sres[w2]))
        testwavs=np.array(testwavs).flatten()
        stds=np.array(stds).flatten()
        if not silence:
            print("Interpolating standard deviations...")
        finterp=interp1d(testwavs,stds*norm,kind='quadratic',fill_value='extrapolate')
        if not silence:
            print("!!!!! FAKING ERROR SPECTRUM !!!!!")
        eflux = finterp(wav)/norm
        w=eflux<np.min(stds)
        if np.sum(w)>0:
            eflux[w] = np.min(stds)

    elif method == "bspl":
        norm=1.00/np.median(flux)
        if lines is None:
            lines=load_lines()
        if not silence:
            print("Fitting an initial model...")
        cspacing=wav/100.0
        cspl,outmask,tck=bspline_iterfit(wav,flux*norm,yerr=None,ivar=None,knots=None,
                                        silence=silence,knotfun=get_contset_placed,
                                        cspacing=cspacing,lines=lines,zguess=zguess)
        if not isinstance(cspl, float):
            fmod=cspl(wav)
        else:
            placed=get_contset_placed(wav,cspacing,lines=lines,zguess=zguess)
            cset=bspline.iterfit(wav,flux*norm,placed=placed)
            fmod=cset[0].value(wav)[0]
        res=flux-fmod/norm
        dwav=np.median(wav[1:]-wav[:-1])
        step=dwav*30.
        testwavs=[];stds=[]
        for thiswav in np.arange(np.min(wav),np.max(wav)-step*0.5,step):
            w=(wav>=thiswav)&(wav<thiswav+step)
            nw=np.sum(w)
            if nw<10:
                continue
            sres=res[w]
            sres=np.sort(sres)
            std=np.std(sres[int(0.1*(nw-1)):int(0.9*(nw-1))])
            medres=np.median(sres)
            w2=np.abs(sres-medres)<(3.0*std)
            testwavs.append(np.mean(wav[w]))
            stds.append(np.std(sres[w2]))
        testwavs=np.array(testwavs).flatten()
        stds=np.array(stds).flatten()
        wavrange=np.max(wav)-np.min(wav)
        kspace=10.0*step
        wavrange=np.arange(testwavs.min(),testwavs.max()+kspace,kspace)
        wavrange=np.delete(wavrange,[0,-1])
        knts=[wav for wav in wavrange]
        if not silence:
            print("Fitting for standard deviations...")
        sspl,outmask,tck=bspline_iterfit(testwavs,stds*norm,knots=None,silence=silence,
                                         knotfun=default_knotfun,kspace=kspace)
        if not silence:
            print("!!!!! FAKING ERROR SPECTRUM !!!!!")
        eflux=sspl(wav)/norm
        w=eflux<np.min(stds)
        nw=np.size(w)
        if nw>0:
            eflux[w]=np.min(stds)

    else:
        print("ESTIMATE_EFLUX: error estimation method not recognized...")
        print("Defaulting to equal weights")
        return np.ones(len(wav))

    return eflux

def convert_double_spec(fname,ecol=2,isest=False):
    """
    FUNCTION: CONVERT_DOUBLE_SPEC
    ****************************************************************************
    Converts a spectrum file to change double exponential format (e.g. 1d-04) to
    one Python can handle
    ****************************************************************************
    INPUT:
        fname   :   Filename for the spectrum | character
    OPTIONAL INPUT:
        ecol    :   Column corresponding to error spectrum | integer, default 2
        isest   :   Set True if error spectrum is estimated | bool, default False
    OUTPUT:
        spec    :   Read-in spectrum | real
    ****************************************************************************
    NOTES:
        -I wrote this because Python doesn't like exponentials of form `XD-XX`
    ****************************************************************************
    Written by Michael Baer
    """
    strspec = np.loadtxt(fname, comments='#', dtype=str).T
    if len(strspec) < 2:
        print('READ_SPEC ERROR: The data file is invalid')
        print('Check README.md for more details')
        return
    dl = []
    for ll in range(len(strspec)):
        for i, dpt in enumerate(strspec[ll]):
            dpt = dpt.replace('D','e')
            strspec[ll,i] = dpt
        dl.append(strspec[ll].astype('float'))
    spec = []
    for i in range(np.size(dl[0])):
        if len(strspec) == 2 or isest:
            spec.append(np.array((dl[0][i],dl[1][i])))
        else:
            spec.append(np.array((dl[0][i],dl[1][i],dl[ecol][i])))
    spec = np.array(spec).T
    return spec

def read_spec(fname, lines=None, zguess=0., isest=False, ecol=2, **kwargs):
    """
    FUNCTION: READ_SPEC
    ****************************************************************************
    Reads a user-provided ASCII-formatted data file for a given spectrum. If
    only 2 columns exist in the file, only wavelength and flux are assumed, with
    an error spectrum generated with a B-Spline iterative fit by default.
    ****************************************************************************
    INPUT:
        fname   :   File name of the given supernova spectrum | character
    OPTIONAL INPUT:
        lines   :   Line information | numpy structure, default None
            lines['ION']       --> Ion name | character
            lines['FIT']       --> True if line is ignored in a fit | bool
            lines['LINE']      --> Wavelength of line in Angstroms | real
            lines['EW']        --> Equivalent width of the line (-1 if
                                       absorption, +1 if emission) | real
            lines['TIEDTO']    --> Something | real
            lines['TIEDRATIO'] --> Something else | real
            lines['ISABS']     --> True if line is absorption | bool
            lines['MASKRAD']   --> Width for line masking in Angstroms | real
        zguess  :   Supernova redshift | real, default 0.0
        isest   :   Even if a third column exists, estimate the error | bool,
                        default False
    OUTPUT:
        outspec :   Structured array containing the loaded spectrum | real

    NOTES:
        If `isest` is True, `method` should be provided as a keyword argument.
            Check ESTIMATE_EFLUX for more details
    ****************************************************************************
    Written by Michael Baer
    """
    # Load the file
    if not os.path.isfile(fname):
        print("READ_SPEC ERROR: file not found")
        return 0., 0.
    try:
        spec = np.loadtxt(fname, comments='#').T
    except ValueError:
        spec = convert_double_spec(fname,ecol=ecol,isest=isest)
        if spec is None:
            return 0., 0.
    #spec = np.loadtxt(fname, comments='#').T
    if len(spec) < 2:
        print('READ_SPEC ERROR: The data file is invalid')
        print('Check README.md for more details')
        return 0., 0.
    elif len(spec) < 3 and not isest:
        print('No flux errors found in data file... estimating!')
        wav = spec[0]
        flux = spec[1]
        wok = np.isfinite(flux)
        wav=wav[wok];flux=flux[wok]
        # If no error spectrum is provided, estimate
        isest = True
        e_flux = estimate_eflux(wav, flux, lines, zguess, **kwargs)
        if not np.any(np.isfinite(e_flux)):
            print("whoops!")
            return 0., 0.
    else:
        if len(spec) == 3:
            ecol = 2
        else:
            ecol = 3
        wav = spec[0]
        flux = spec[1]
        if isest:
            wok = np.isfinite(flux)
            wav=wav[wok];flux=flux[wok]
            e_flux = estimate_eflux(wav, flux, lines, zguess, **kwargs)
        else:
            e_flux = spec[ecol]
            if np.sum(e_flux > 0.) == 0:
                isest = True
                wok = np.isfinite(flux)
                wav=wav[wok];flux=flux[wok]
                e_flux = estimate_eflux(wav, flux, lines, zguess, **kwargs)
    # Clean up
    w = (np.isfinite(flux) & np.isfinite(e_flux))
    wav = wav[w]
    flux = flux[w]
    e_flux = e_flux[w]
    # Save as structured array
    outspec = []
    for i in range(np.size(wav)):
        outspec.append(np.array((wav[i], flux[i], e_flux[i]), dtype=[('wav', float), ('flux', float), ('eflux', float)]))
    outspec = np.array(outspec)
    return outspec, isest

def bin_spec_weighted(spec,lowav,hiwav,dwav,sigma=None,samewav=True,sigwav=None):
    """
    FUNCTION: BIN_SPEC_WEIGHTED
    ****************************************************************************
    Bins a given spectrum, using weights if provided
    ****************************************************************************
    INPUT:
        spec    :   Supernova spectrum | real
            spec['wav']     --> Wavelengths
            spec['flux']    --> Flux densities
        lowav   :   Minimum wavelength | real
        hiwav   :   Maximum wavelength | real
        dwav    :   Bin width | real
    OPTIONAL INPUT:
        sigma   :   User-provided weights; if not provided, defaults to equal
                        weights | real, default None
        samewav :   Set True if wavelengths of spectrum and weights are the
                        same | bool, default True
        sigwav  :   Wavelengths of weight spectrum | real, default None
    OUTPUT:
        outspec :   Structured array containing the binned spectrum | real
    ****************************************************************************
    Written by Michael Baer
		-adapted from Superfit in IDL by Andy Howell, written to utilize numpy
			structured arrays
    """
    if sigma is not None:
        if not samewav:
            if sigwav is None:
                raise ValueError("Input weights need wavelengths")
        else:
            sigwav=np.copy(spec['wav'])
    else:
        sigma=np.ones(len(spec['wav']))
        sigwav=np.copy(spec['wav'])
    if dwav==0.:
        nwav=len(spec['wav'])
        outspec=np.zeros(nwav,dtype=[('wav',float),('flux',float),('weight',float)])
        wh=(spec['wav']>=lowav)&(spec['wav']>=hiwav)
        outspec['wav']=spec['wav'][wh]
        outspec['flux']=spec['flux'][wh]
        if not samewav:
            raise ValueError("Weight and flux must be same size to rescale at same resolution")
        outspec['weight']=sigma
    elif dwav>0.:
        nwav=int((hiwav-lowav)/dwav)+1
        outspec=np.zeros(nwav,dtype=[('wav',float),('flux',float),('weight',float)])
        outspec['wav']=np.arange(nwav)*dwav+lowav+dwav/2.0
        interpwav=np.append(spec['wav'],outspec['wav'])
        interpwav=np.sort(interpwav)
        interpwav=np.unique(interpwav)
        finterp=interp1d(spec['wav'],spec['flux'],bounds_error=False)
        interpflux=finterp(interpwav)
        wok=np.isfinite(interpflux)
        ginterp=interp1d(sigwav,sigma,bounds_error=False)
        interpsig=ginterp(interpwav)
        wok=np.isfinite(interpsig)

        for i in range(0,nwav):
            w=(interpwav>=outspec['wav'][i]-dwav/2.0)&(interpwav<=outspec['wav'][i]+dwav/2.0)
            outspec['flux'][i]=np.nansum(interpflux[w]/interpsig[w]**2)/np.nansum(1./interpsig[w]**2)
            outspec['weight'][i]=np.nanmean(interpsig[w])
    return outspec

def cleanup_spec(spec,niter,grow,sigma=None,method='incl',nsig=3):
    """
    FUNCTION: CLEANUP_SPEC
    ****************************************************************************
    Performs sigma-clipping on a spectrum
    ****************************************************************************
    INPUT:
        spec    :   Supernova spectrum | real
            spec['wav']     --> Wavelengths
            spec['flux']    --> Flux densities
        niter   :   Number of iterations | integer
        grow    :   Grow radius; if greater than 0, replaces values on [grow]
                        pixels around bad pixels | integer
    OPTIONAL INPUT:
        sigma   :   User-provided weights; if not provided, defaults to equal
                        weights | real, default None
        method  :   Method of calculating sigma | character, default `incl`
            `incl`  --> Included (sigma must not be None)
            `calc`  --> Calculate as running standard deviation
        nsig    :   Multiple of sigma past which to clip | real, default 3
    OUTPUT:
        outspec :   Structured array containing the cleaned spectrum | real
    ****************************************************************************
    Written by Michael Baer
		-adapted from Superfit in IDL by Andy Howell, written to utilize numpy
			structured arrays
    """
    nwav=len(spec['wav'])
    fixedflux=np.copy(spec['flux'])
    if method=='incl':
        if sigma is None:
            raise ValueError("incl method needs provided sigma")
    for j in range(niter):
        md=median_filter(fixedflux,size=20)
        if method=='stdev':
            sigma=np.zeros(nwav)
            for i in range(10,nwav-11):
                sigma[i]=np.std(fixedflux[i-10:i+10])
            sigma[:9]=sigma[10]
            sigma[-10:]=sigma[-11]
        w=(fixedflux>md+nsig*sigma)|(fixedflux<md-nsig*sigma)
        for i in range(grow+1):
            fixedflux[w+i]=md[w+i]
            fixedflux[w-i]=md[w-i]
    outspec=np.copy(spec)
    outspec['flux']=fixedflux
    return outspec
