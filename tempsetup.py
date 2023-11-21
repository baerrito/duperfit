"""
    SCRIPT: TEMPSETUP
    *************************************
    Sets up and saves template pickle files

    May be used to implement custom templates; however, I recommend giving a new
    `allsne` dictionary a new name to avoid unexpected bugs
"""

# Import modules
import os
import glob
import pickle
import numpy as np
# Custom packages
from spec_handling import read_spec

try:
    dfdir=os.environ['DFDIR']
except KeyError:
    raise Exception("$DFDIR not recognized")

def setup_template_dictionary(thetype,tempdir):
    tfiles=glob.glob(tempdir+f'{thetype}/*.dat')
    temdict={}
    temdict['tfiles']=[]
    temdict['tspecs']=[]
    for tfile in tfiles:
        tspec,isest=read_spec(tfile)
        temdict['tfiles'].append(tfile)
        temdict['tspecs'].append(tspec)
    with open(f"./picklejar/{thetype}.pickle","wb") as f:
        pickle.dump(temdict, f)
    f.close()
    return

def setup_allsne_dictionary(tempdir,name="allsne"):
    tfiles=glob.glob(tempdir + "*/*.dat")
    temdict={}
    temdict['tfiles']=[]
    temdict['tspecs']=[]
    for tfile in tfiles:
        tspec,isest = read_spec(tfile)
        temdict['tfiles'].append(tfile)
        temdict['tspecs'].append(tspec)
    with open(f"./picklejar/{name}.pickle","wb") as f:
        pickle.dump(temdict, f)
    f.close()

def setup_gal_dictionary(galdir):
    gfiles=glob.glob(galdir + "*")
    galdict={}
    galdict['gfiles']=[]
    galdict['gspecs']=[]
    for gfile in gfiles:
        gspec,isest = read_spec(gfile)
        galdict['gfiles'].append(gfile)
        galdict['gspecs'].append(gspec)
    with open(f"./picklejar/gal.pickle","wb") as f:
        pickle.dump(galdict, f)
    f.close()

if __name__ == '__main__':

    tempdir='./tempbank/sne/'
    galdir='./tempbank/gal/'
    thedirs=glob.glob(tempdir+'*')
    thetypes=[]
    for direc in thedirs:
        typ=direc[direc.rfind("/")+1:]
        #if "SLSN" not in typ:
        #    typ=typ.split("-")[0]
        thetypes.append(typ)
    typetrunc=[*set(thetypes)]

    setup_allsne_dictionary(tempdir)
    for thetype in typetrunc:
        setup_template_dictionary(thetype,tempdir)
    setup_gal_dictionary(galdir)
