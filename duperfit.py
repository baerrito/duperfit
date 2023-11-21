"""
    MODULE: DUPERFIT
    *************************************
    The main module of the code, contains all core functions

	I tried to document each function as thoroughly as possible, but feedback is
	always appreciated
"""

# Import modules
import os
import glob
import time
import pickle
import numpy as np
from astropy.table import Table
import astropy.io.ascii as tasc
from scipy.optimize import curve_fit, minimize
import extinction
import warnings
# Custom modules
from spec_handling import *
# Plotting
import matplotlib.pyplot as plt
plt.rc('axes', labelsize=14)
plt.rc('axes', labelweight='bold')
plt.rc('figure', titlesize=16)
plt.rc('figure', titleweight='bold')
plt.rc('font', family='sans-serif')
plt.rcParams['errorbar.capsize'] = 3
opts = {'mec':'k', 'mew': 0.5, 'lw': 1}
plt.rcParams['figure.figsize'] = (20, 10)
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True

warnings.filterwarnings("ignore")

# Functions
def sf_model(X,*params):
	"""
	FUNCTION: SF_MODEL
	****************************************************************************
	Model function
	****************************************************************************
	INPUT:
		X			:	Sequence of dependent variable arrays | real
			X[0]	-->		Template flux
			X[1]	-->		Galaxy flux
			X[2]	-->		Reddening law divided by A_V
		params		:	Sequence of fitted parameters | real
			params[0]	-->	Template scaling
			params[1]	-->	Galaxy scaling
			params[2]	-->	A_V
	OUTPUT:
		sflux		:	Model supernova flux | real
	****************************************************************************
	Written by Michael Baer
		-adapted from Superfit in IDL by Andy Howell
	"""
	tflux,gflux,redlaw = X
	cc,ff,aa = params
	sflux = cc*tflux*10**(aa*redlaw) + ff*gflux
	return sflux

# Functions
def sf_model_fixg(X,ff,*params):
	"""
	FUNCTION: SF_MODEL_FIXG
	****************************************************************************
	Model function with fixed galaxy flux
	Written because scipy.optimize doesn't play well with parameters with equal
		lower/upper bounds
	****************************************************************************
	INPUT:
		X			:	Sequence of dependent variable arrays | real
			X[0]	-->		Template flux
			X[1]	-->		Galaxy flux
			X[2]	-->		Reddening law divided by A_V
		ff			:	Galaxy scaling (fixed) | real
		params		:	Sequence of fitted parameters | real
			params[0]	-->	Template scaling
			params[1]	-->	A_V
	OUTPUT:
		sflux		:	Model supernova flux | real
	****************************************************************************
	Written by Michael Baer
		-adapted from Superfit in IDL by Andy Howell
	"""
	tflux,gflux,redlaw = X
	cc,aa = params
	sflux = cc*tflux*10**(aa*redlaw) + ff*gflux
	return sflux

def optimize_func(params,X,oflux,weight):
	"""
	FUNCTION: OPTIMIZE_FUNC
	****************************************************************************
	Target function for optimization
	****************************************************************************
	INPUT:
		params		:	Sequence of fitted parameters | real
		X			:	Sequence of dependent variable arrays | real
			X[0]	-->		Template flux
			X[1]	-->		Galaxy flux
			X[2]	-->		Reddening law divided by A_V
		oflux		:	Object flux | real
		weight		:	Weight array | real
	OUTPUT:
		wgresid		:	Weighted sum of square residuals | real
	****************************************************************************
	Written by Michael Baer
		-adapted from Superfit in IDL by Andy Howell
	"""
	omod=sf_model(X,*params)
	wgresid=np.nansum( (oflux-omod)**2*weight )
	return wgresid

def optimize_func_fixg(params,X,ff,oflux,weight):
	"""
	FUNCTION: OPTIMIZE_FUNC_FIXG
	****************************************************************************
	Target function for optimization, with fixed galaxy flux
	Written because scipy.optimize doesn't play well with parameters with equal
		lower/upper bounds
	****************************************************************************
	INPUT:
		params		:	Sequence of fitted parameters | real
		X			:	Sequence of dependent variable arrays | real
			X[0]	-->		Template flux
			X[1]	-->		Galaxy flux
			X[2]	-->		Reddening law divided by A_V
		ff			:	Galaxy scaling (fixed) | real
		oflux		:	Object flux | real
		weight		:	Weight array | real
	OUTPUT:
		wgresid		:	Weighted sum of square residuals | real
	****************************************************************************
	Written by Michael Baer
		-adapted from Superfit in IDL by Andy Howell
	"""
	omod=sf_model_fixg(X,ff,*params)
	wgresid=np.nansum( (oflux-omod)**2*weight )
	return wgresid

def core_fit(z,ospec,tspec,gspec,minw,maxw,resolution,Avmin,Avmax,exactz=False,\
			 sfractol=0.7,Rv=3.1,tscale=3.0,gscale=3.0,silence=False,
			 optimizer="L-BFGS-B"):
	"""
	FUNCTION: CORE_FIT
	****************************************************************************
	Core fitting function for Duperfit, fitting for the model:
			o(lam) = C*s(lam,z)*10^(-0.4*A_V*r(lam)) + D*g(lam,z)
	****************************************************************************
	INPUT:
		z			:	Redshift | real
		ospec		:	Object spectrum | structured array of reals
			ospec['wav']	-->		Wavelengths of the spectrum
			ospec['flux']	-->		Flux densities of the spectrum
			ospec['weight']	-->		Fitting weights of the spectrum (1/sigma^2)
		tspec		:	SN template spectrum | structured array of reals
			tspec['wav']	-->		Wavelengths of the spectrum
			tspec['flux']	-->		Flux densities of the spectrum
		gspec		:	Galaxy template spectrum | structured array of reals,
							same structure as tspec
		minw		:	Minimum wavelength for fit | real
		maxw		:	Maximum wavelength for fit | real
		resolution	:	Wavelength resolution for fit | real
		Avmin 		:	Lower bound on V-band extinction | real
		Avmax 		:	Upper bound on V-band extinction | real
	OPTIONAL INPUT:
		exactz		:	Set True if z is fixed | bool, default False
		sfractol	:	Lower bound on fractional wavelength coverage | real,
							default 0.7
		Rv			:	Cardelli extinction law parameter | real, default 3.1
		tscale		:	Upper bound on SN scaling C | real, default 3.0
		gscale 		:	Upper bound on Gal scaling D | real, default 3.0
		silence		:	Set True to silence messages | bool, default False
		optimizer	:	Choice of optimizer (L-BFGS-B with minimize or TRF
							with curve_fit) | character, default `L-BFGS-B`
	OUTPUT:
		S 			:	Optimization value | real
		popt		:	Optimal values of C, D, A_V | array of reals
		sfrac 		:	Fractional wavelength coverage by SN template | real
	****************************************************************************
	NOTES:
		-`ospec` should be binned such that minw<=spec['wav']<=maxw on input
		-might add option to use either minimize with L-BFGS-B or curve_fit with
			Trust Region Reflective or Dog-box... still deciding
	****************************************************************************
	Written by Michael Baer
		-adapted from Superfit in IDL by Andy Howell, written to utilize numpy
			structured arrays
	"""
	objspec = np.copy(ospec)
	sflag=0
	twavred=tspec['wav']*(1+z)
	if minw<twavred[0]:
		sminw=twavred[0]
		sflag=1
	else:
		sminw=minw
	if maxw>twavred[-1]:
		smaxw=twavred[-1]
		sflag=1
	else:
		smaxw=maxw
	tspecz=np.copy(tspec)
	tspecz['flux']/=np.median(tspecz['flux'])
	tspecz['wav']=twavred
	gspecz=np.copy(gspec)
	gspecz['flux']/=np.median(gspecz['flux'])
	gspecz['wav']=gspec['wav']*(1+z)
	sfrac=(smaxw-sminw)/(maxw-minw)
	if sfrac<sfractol:
		if not silence:
			print("Wavelength coverage insufficient... not fitting")
		popt=np.zeros(3)
		pcov=np.zeros((3,3))
		return np.inf, popt, sfrac, np.inf
	glspec = bin_spec_weighted(gspecz,sminw,smaxw,resolution)
	tmspec = bin_spec_weighted(tspecz,sminw,smaxw,resolution)
	if sflag==1:
		sobjspec=bin_spec_weighted(objspec,sminw,smaxw,resolution,objspec['weight'])
		noredw=sobjspec['wav']/(1+z)
		objflux=sobjspec['flux']
		objweight=sobjspec['weight']
	else:
		noredw=objspec['wav']/(1+z)
		objflux=objspec['flux']
		objweight=objspec['weight']

	# compute r(lam)
	redlaw = extinction.ccm89(noredw,1.,Rv)

	# define the X-tuple for the model function
	Xfit = (tmspec['flux'],glspec['flux'],redlaw)
	guess = [0.6,0.6,0.0]
	aamin=-0.4*Avmax
	aamax=-0.4*Avmin

	if optimizer=='TRF':
		# Method 1 (slower, can give uncertainties on pars, though errors are poor)
		# May eliminate this in future version
		bounds = ([0.01,0.0,aamin],[tscale,gscale,aamax])
		try:
			popt, pcov = curve_fit(sf_model,Xfit,objflux,p0=guess,bounds=bounds,\
						sigma=1./np.sqrt(objweight))
		except RuntimeError:
			if not silence:
				print("Optimal parameters not found by curve_fit")
			popt=np.zeros(3)
			pcov=np.zeros((3,3))
			return np.inf, popt, sfrac, np.inf
		omod = sf_model(Xfit,*popt)
		S = np.nansum((objspec['flux']-omod)**2*objspec['weight'])/sfrac

	elif optimizer=='L-BFGS-B':
		# Method 2 (faster, equally good fits)
		if gscale==0.0:
			guess = [0.6,0.0]
			bounds=[(0.01,tscale),(aamin,aamax)]
			res = minimize(optimize_func_fixg,guess,args=(Xfit,0.0,objflux,objweight,),bounds=bounds,method=optimizer)
		else:
			bounds=[(0.01,tscale),(0.0,gscale),(aamin,aamax)]
			res = minimize(optimize_func,guess,args=(Xfit,objflux,objweight,),bounds=bounds,method=optimizer)
		if not res.success:
			if not silence:
				print("Optimal parameters not found by minimize")
			popt=np.zeros(3)
			return np.inf, popt, sfrac, np.inf
		popt = res.x
		#ddof=np.size(np.logical_not())
		if gscale==0.0:
			popt=[popt[0],0.0,popt[1]]
		modflux=sf_model(Xfit,*popt)
		ndat=np.sum(np.logical_not(np.isnan((objflux-modflux)**2*objweight)))
		if exactz:
			ddof=ndat-3
		else:
			ddof=ndat-4
		S = res.fun/sfrac

	else:
		raise NameError(f"optimizer='{optimizer}' not recognized")

	# gfrac calculation; really not too expensive, just set infinite if no
	#   optimization was yielded
	if np.isfinite(S):
		galsum=0.;obssum=0.
		gterm=gspecz['flux']*popt[1]
		# Approximate V-band wavelength range
		gfracwavstart=4000*(1+z)
		gfracwavend=5000*(1+z)

		# Sum the galaxy flux
		wx=(gspecz['wav']>=gfracwavstart)&(gspecz['wav']<=gfracwavend)
		if np.sum(wx)>0:
			galsum=np.sum(gterm[wx])
		else:
			galsum=np.sum(gterm)

		# Sum the SN flux
		if sflag==1:
			wy=(sobjspec['wav']>=gfracwavstart)&(sobjspec['wav']<=gfracwavend)
		else:
			wy=(objspec['wav']>=gfracwavstart)&(objspec['wav']<=gfracwavend)
		if np.sum(wy)>0:
			obssum=np.sum(objflux[wy])
		else:
			obssum=np.sum(objflux)

		gfrac=galsum/obssum

	else:

		gfrac=np.inf

	return S, popt, sfrac, gfrac

def duperfit(zmin,zmax,zstep,ospec,objname,temdict,galdict,Avmin=-2.0,Avmax=2.0,
			 wgdir=None,wsrc='uw',sclip=True,sigsrc='calc',niter=5,grow=0,
			 nsig=2.7,exactz=False,beginw=2000,endw=10000,resolution=20.,
			 sfractol=0.7,Rv=3.1,tscale=3.0,gscale=3.0,silence=False,save=False,
			 output_path="",outfile="",optimizer="L-BFGS-B"):
	"""
	FUNCTION: DUPERFIT
	****************************************************************************
	Performs Superfit-like optimization
	****************************************************************************
	INPUT:
		zmin		:	Minimum evaluating redshift | real
		zmax		:	Maximum evaluating redshift | real
		zstep		:	Grid size for evaluating redshift | real
		ospec		:	Object spectrum | structured array of reals
			ospec['wav']	-->		Wavelengths of the spectrum
			ospec['flux']	-->		Flux densities of the spectrum
			ospec['eflux']	-->		Flux uncertainties [optional]
		objname		:	Name of the object | character
		temdict		:	SN template dictionary | dictionary
			temdict['tfiles']
				-->	Filenames for SN templates | character
			temdict['tspecs']
				-->	SN template spectra | structured arrays of reals
					tspec['wav']	-->		Wavelengths of the spectrum
					tspec['flux']	-->		Flux densities of the spectrum
		galdict		:	Gal template dictionary | dictionary
			galdict['gfiles']
				-->	Filenames for Gal templates | character
			galdict['gspecs']
				-->	Gal template spectra | structured arrays of reals
					gspec['wav']	-->		Wavelengths of the spectrum
					gspec['flux']	-->		Flux densities of the spectrum
	OPTIONAL INPUT:
		Avmin 		:	Lower bound on V-band extinction | real, default -2.0
		Avmax 		:	Upper bound on V-band extinction | real, default 2.0
		wgdir		:	Directory containing weight file | character, default None
		wsrc		:	Source of weights | character, default `uw`
			`uw`		-->	Unweighted (1.0 at all points)
			`tell`		-->	Telluric deweighted
			`incl`		-->	Included
		sclip		:	Set True if ospec should be sigma clipped | bool,
							default True
		sigsrc		:	Source for values for sigma clipping | character,
							default `calc`
			`calc`		-->	Calculate standard deviations
			`weight`	--> Use weights
		niter		:	Number of iterations for sigma clipping | integer,
							default 5
		grow		:	Grow parameter for sigma clipping | integer, default 0
		nsig 		:	Number of sigmas to sigma clip | real, default 2.7
		exactz		:	Set True if z is fixed | bool, default False
		beginw		:	Minimum wavelength for fit | real, default 2000
		endw		:	Maximum wavelength for fit | real, default 10000
		resolution	:	Wavelength resolution for fit | real, default 10.0
		sfractol	:	Lower bound on fractional wavelength coverage | real,
							default 0.7
		Rv			:	Cardelli extinction law parameter | real, default 3.1
		tscale		:	Upper bound on SN scaling C | real, default 3.0
		gscale 		:	Upper bound on Gal scaling D | real, default 3.0
		silence		:	Set True to silence messages | bool, default False
		save		:	Set True to save output	| bool, default False
		output_path	:	Directory for output file | character, default ``
							(current directory)
		optimizer	:	Choice of optimizer (L-BFGS-B with minimize or TRF
							with curve_fit) | character, default `L-BFGS-B`
	OUTPUT:
		output 		:	Astropy table of outputs | table
			output['SN']
				-->	File paths for SN templates | character
			output['S']
				-->	S value | real
			output['z']
				--> Redshift | real
			output['gal']
				-->	Galaxy types | character
			output['cSN']
				--> C | real
			output['cgal']
				--> D | real
			output['Av']
				--> A_V | real
			output['sfrac']
				--> Fractional wavelength coverage | real
	****************************************************************************
	NOTES:
		-`ospec['eflux']` is required if `wsrc='incl'`.
		-`output` is sorted from least to greatest chi^2
		-`sigsrc='weight'` is not recommended; still needs work
		-`wsrc='calc'` may eventually be implemented; if so, recommended for
			consistency (either B-spline iterative fit or SG filter)
	****************************************************************************
    Written by Michael Baer
		-adapted from Superfit in IDL by Andy Howell, written to utilize numpy
			structured arrays
	"""

	# Check wavelength range
	minw=np.min(ospec['wav'])
	maxw=np.max(ospec['wav'])

	# Raise error if spectrum isn't viable
	if maxw<beginw:
		raise ValueError("Duperfit failed; input spectrum below minimum wavelength coverage")
	if minw>endw:
		raise ValueError("Duperfit failed; input spectrum above maximum wavelength coverage")

	# Trim min/max wavelengths to the spectrum, if necessary
	if beginw>=minw:
		minw=beginw
	if endw<=maxw:
		maxw=endw

	# Get weights
	if wsrc=='incl':
		if 'eflux' not in ospec.dtype.names:
			raise ValueError("Duperfit failed; input spectrum lacking error fluxes when wsrc=`incl`")
		weight=ospec['eflux']

	elif wsrc=='tell':
		#if wgdir is None:
		#	raise NameError('wgdir not provided')
		wgspec=np.loadtxt('./tempbank/no77.weight').T
		weight=wgspec[1]
		wgwav=wgspec[0]
	elif wsrc=='uw':
		weight=np.ones(len(ospec['wav']))

	# Clean up the spectrum, if called for
	if sclip:
		if sigsrc=='weight':
			ospec=cleanup_spec(ospec,niter,grow,weight,nsig=nsig)
		elif sigsrc=='calc':
			ospec=cleanup_spec(ospec,niter,grow,method='stdev',nsig=nsig)

	# Normalize stuff
	weight=weight/np.median(weight)
	ospec['flux']=ospec['flux']/np.median(ospec['flux'])

	# Enforce 1/sigma^2
	weight = 1/weight**2

	# Enforce wavelength bounds are doable by object spectrum
	if minw<ospec['wav'][0]:
		minw=ospec['wav'][0]
	if maxw>ospec['wav'][-1]:
		maxw=ospec['wav'][-1]

	# Bin the object spectrum to wavelength bounds
	if wsrc=='tell':
		bospec = bin_spec_weighted(ospec,minw,maxw,resolution,weight,samewav=False,sigwav=wgwav)
	else:
		bospec = bin_spec_weighted(ospec,minw,maxw,resolution,weight)

	# Get number of redshift steps
	if exactz:
    		nzs=1
	else:
    		if zstep<=0:
    			print("zstep not set, defaulting to 0.01")
    			zstep=0.01
    		nzs = int((zmax-zmin)/zstep)+1

	# Get file names and template spectra
	# This is very important, I promise
	tfiles=temdict['tfiles']
	tspecs=temdict['tspecs']
	gfiles=galdict['gfiles']
	gspecs=galdict['gspecs']

	# Initialize some arrays cuz that's cool
	#S=np.zeros((len(tfiles),len(gfiles),nzs))
	tbS=np.zeros(len(tfiles))
	tbG=np.zeros(len(tfiles),dtype='<U3')
	tbz=np.zeros(len(tfiles))
	tbcc=np.zeros(len(tfiles))
	tbff=np.zeros(len(tfiles))
	tbAv=np.zeros(len(tfiles))
	#tbae=np.zeros(len(tfiles))
	#tbbe=np.zeros(len(tfiles)) ## Pay no attention to the parameter errors behind the curtain
	#tbAve=np.zeros(len(tfiles))
	tbsf=np.zeros(len(tfiles))
	tbgf=np.zeros(len(tfiles))

	# Initialize some values so Python doesn't scream at me
	#bestzfort=0.
	#bestG=0
	#bestGa=0.
	#bestGb=0.
	#bestGAv=0.
	##bGaerr=0.
	##bGberr=0.
	##bGAverr=0.
	#bGsfrac=0.
	#bGgfrac=0.

    # Initializing arrays for minimization
    S=np.zeros((len(gfiles),nzs))
    Gal=np.zeros((len(gfiles),nzs),dtype='<U3')
    redshift=np.zeros((len(gfiles),nzs))
    cc=np.zeros((len(gfiles),nzs))
    ff=np.zeros((len(gfiles),nzs))
    Av=np.zeros((len(gfiles),nzs))
    #ccerr=np.zeros((len(gfiles),nzs))
    #fferr=np.zeros((len(gfiles),nzs))
    #Averr=np.zeros((len(gfiles),nzs))
    sfrac=np.zeros((len(gfiles),nzs))
    gfrac=np.zeros((len(gfiles),nzs))
    S[:,:]=np.inf
    gfrac[:,:]=np.inf

	# Initialize a template list... this is also very important
	thetemps=[]

	# Tell the users what's happening, it's probably a good idea
	print(f"Running for {objname} at resolution = {resolution:.0f} Angstrom")
	# Load bars are pretty cool, we should include one of those
	if silence:
		print(f'Computing: [{"."*60}] {0:.1f}%',end='\r',flush=True)

	# SN loop
	for i,(tfile,tspec) in enumerate(zip(tfiles,tspecs)):
		parts=tfile.split("/")
		tname=parts[-2]+"/"+parts[-1]
		thetemps.append(tname)

		# Galaxy loop
		for j,(gfile,gspec) in enumerate(zip(gfiles,gspecs)):
			gname=gfile[gfile.rfind("/")+1:]
			#if j==0:
			#	bestGS=np.inf
			#	bestGa=0.
			#	bestGb=0.
			#	bestGAv=0.
			#	bGgfrac=np.inf

			# Redshift loop
			for k in range(nzs):
				if k==0:
					bestzS=np.inf
					bzgfrac=np.inf

				# Set z
				z=zmin+k*zstep
                redshift[j,k]=z

				# Run the core function
				S[j,k],popt,sfrac[j,k],gfrac[j,k]=core_fit(z,bospec,tspec,\
                                  gspec,minw,maxw,resolution,Avmin,Avmax,\
                                  exactz=exactz,sfractol=sfractol,Rv=Rv,\
                                  tscale=tscale,gscale=gscale,silence=silence)
                cc[j,k]=popt[0]
                ff[j,k]=popt[1]
                Av[j,k]=-2.5*popt[2]

			# End of the redshift loop

		#End of the galaxy loop

        # Determine loci with minimum S
        wmin=np.unravel_index(np.argmin(S),S.shape)
		tbS[i]=S[wmin]
		tbG[i]=Gal[wmin]
		tbz[i]=redshift[wmin]
		tbcc[i]=cc[wmin]
		tbff[i]=ff[wmin]
		tbAv[i]=Av[wmin]
		#tbae[i]=bGaerr
		#tbbe[i]=bGberr
		#tbAve[i]=bGAverr
		tbsf[i]=sfrac[wmin]
		tbgf[i]=gfrac[wmin]

		# Updating the load bar
		percent=(i+1)*100.0/len(tfiles)
		x=int(60*percent/100)
		bleep='#'*x;bloop="."*(60-x) # I have the worst variable names sometimes
		if silence:
			print(f'Computing: [{bleep}{bloop}] {percent:.1f}%',end='\r',flush=True)

	if silence:
		print('\n',flush=True) # Line break after load bar so terminal output isn't uggo
	# End of the loop

    # Deallocate unneeded arrays
    del S, Gal, redshift, cc, ff, Av, sfrac, gfrac

	# We're done! ... Or are we?
	thetemps=np.array(thetemps) # *Vsauce music*

	# Sorting by increasing S (better to worse fits)
	wsort=np.argsort(tbS)

	# This output table never happened
	#inds=np.arange(np.size(wsort))
	#output = table.Table((inds,thetemps[wsort],tbG[wsort],tbS[wsort],tbz[wsort],\
	#                      tba[wsort],tbae[wsort],tbb[wsort],tbbe[wsort],tbAv[wsort],\
	#                      tbAve[wsort],tbsf[wsort]),
	#                     names=('I','SN','gal','S','z','cSN','cSNerr','cgal','cgalerr',\
	#                            'Av','Averr','sfrac'),
	#                     dtype=('int','str','str','float','float','float','float','float','float',\
	#                            'float','float','float'))

	# This output table actually happened
	output = Table((thetemps[wsort],tbS[wsort],tbz[wsort],tbG[wsort],\
			       tbAv[wsort],tba[wsort],tbb[wsort],tbsf[wsort],tbgf[wsort]),\
			      names=('SN','S','z','gal','Av','cSN','cgal','sfrac','gfrac'),\
			      dtype=('str','float','float','str','float','float','float',\
				  		 'float','float'))

	# Updating the formats to be reasonably readable by humans
	output['S'].info.format='.3f'
	output['z'].info.format='.4f'
	output['Av'].info.format='.4f'
	output['cSN'].info.format='.4f'
	output['cgal'].info.format='.4f'
	output['sfrac'].info.format='.3f'
	output['gfrac'].info.format='.3f'

	# Save a file it's really cool
	if save:
		if outfile=="":
			outfile=f"sf_result_{objname}.dfo"
		fname=output_path+outfile
		print(f"Saving output to {fname}") # Yeah it's probably good if someone knows that a file is being saved
		tasc.write(output,fname,format='fixed_width',overwrite=True)

	# Now we're done, here's your table
	return output
