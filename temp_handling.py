import os
import glob
import pickle
import numpy as np

try:
    dfdir=os.environ['DFDIR']
except KeyError:
    raise Exception("$DFDIR not recognized")

def get_tempdict(tempchoice,uinpfile=""):
	"""
	FUNCTION: GET_TEMPDICT
	****************************************************************************
	Function for obtaining a dictionary of SN templates
	****************************************************************************
	INPUT:
		tempchoice	:	User choice of template library | character
							`User Input` will expect a filename
							`SNe <=10d` will use all SN templates with phases
								less than or equal to 10 days from maximum
							`SNIa` truncates all normal SN Ia libraries
							All other choices are filenames <tempchoice>.pickle
								in the `picklejar` directory
	OPTIONAL INPUT:
		uinpfile	:	User-provided template file | character
							Template file needs to be .pickle format, following
								temdict structure (see README.md)
	OUTPUT:
		tempdict	:	SN template dictionary | dictionary
			tempdict['tfiles']
				-->	Filenames for SN templates | character
			tempdict['tspecs']
				-->	SN template spectra | structured arrays of reals
					tspec['wav']	-->		Wavelengths of the spectrum
					tspec['flux']	-->		Flux densities of the spectrum
	****************************************************************************
	Written by Michael Baer
	"""
	if tempchoice=="User Input":
		fname=uinpfile
		if not os.path.isfile(uinpfile):
			raise OSError("File not found")
	elif tempchoice=="SNe <=10d" or tempchoice=="SNIa":
		fname=dfdir+"/picklejar/allsne.pickle"
	else:
		fname=dfdir+f"/picklejar/{tempchoice}.pickle"

	if tempchoice=="SNe <=10d":
		with open(fname,"rb") as f:
			allsne=pickle.load(f)
		f.close()
		allsfs=np.array(allsne['tfiles'])
		allsss=np.array(allsne['tspecs'],dtype=object)
		tempdict={}
		tempdict['tfiles']=[]
		tempdict['tspecs']=[]
		for sf,ss in zip(allsfs,allsss):
			parts=sf.split(".")
			ph=parts[-2]
			if ph=="max":
				phase=0.
			elif ph[0]=="m":
				phase=-1.*float(ph[1:])
			elif ph[0]=="p":
				phase=float(ph[1:])
			if phase<=10:
				tempdict['tfiles'].append(sf)
				tempdict['tspecs'].append(ss)

	elif "SNe" in tempchoice:
		parts=tempchoice.split(" ")
		if parts[-1]=="IIn":
			allfiles=[dfdir+"/picklejar/SNIIn.pickle",dfdir+"/picklejar/SLSN-IIn.pickle",
					  dfdir+"/picklejar/SNIa-CSM.pickle"]
		else:
			allfiles=glob.glob(dfdir+f"/picklejar/SN{parts[-1]}-*.pickle")
		tempdict={}
		tempdict['tfiles']=[]
		tempdict['tspecs']=[]
		for fname in allfiles:
			if "n" not in parts[-1]:
				if "CSM" in fname or "n" in fname:
					continue
			with open(fname,"rb") as f:
				ansndict=pickle.load(f)
			f.close()
			tempdict['tfiles']+=ansndict['tfiles']
			tempdict['tspecs']+=ansndict['tspecs']

	else:
		with open(fname,"rb") as f:
			tempdict=pickle.load(f)
		f.close()

	return tempdict

def get_galdict(galpicks,uinpfile=""):
	"""
	FUNCTION: GET_GALDICT
	****************************************************************************
	Function for obtaining a dictionary of galaxy templates
	****************************************************************************
	INPUT:
		galpicks	:	Sequence of user choices of galaxy templates | character
							`User Input` in the sequence will expect a filename
	OPTIONAL INPUT:
		uinpfile	:	User-provided template file | character
							Template file needs to be an ASCII format galaxy
								spectrum (see README.md)
	OUTPUT:
		galdict		:	Galaxy template dictionary | dictionary
			galdict['gfiles']
				-->	Filenames for galaxy templates | character
			galdict['gspecs']
				-->	Galaxy template spectra | structured arrays of reals
					gspec['wav']	-->		Wavelengths of the spectrum
					gspec['flux']	-->		Flux densities of the spectrum
	****************************************************************************
	Written by Michael Baer
	"""
	with open(dfdir+f"/picklejar/gal.pickle","rb") as f:
		allgaldict=pickle.load(f)
	f.close()
	allgfs=np.array(allgaldict['gfiles'])
	allgss=np.array(allgaldict['gspecs'],dtype=object)
	galdict={}
	galdict['gfiles']=[]
	galdict['gspecs']=[]
	if len(galpicks)>0:
		for gp in galpicks:
			galfile=f"./tempbank/gal/{gp}"
			w=allgfs==galfile
			galdict['gfiles'].append(allgfs[w][0])
			galdict['gspecs'].append(allgss[w][0])
	if "User Input" in galpicks and uinpfile!="":
		uinpspec=np.genfromtxt(uinpfile,names='wav,flux')
		galdict['gfiles'].append(uinpfile)
		galdict['gspecs'].append(uinpspec)
	return galdict
