"""
    SCRIPT: TEMP_MIDSCORES
    *************************************
    Computes MID scores for the template library

    NOTE: uses old-version CSV files which should be included with your
    installation; the score pickle files are also included, this is more meant
    for troubleshooting
"""

# Import modules
import os
import sys
import glob
import pickle
import json
import numpy as np
from astropy import table
import astropy.io.ascii as tasc
# Custom modules
from auxiliary import wheretype

def main(sntypes,basetype):

    basedir="./tempbank/score_eval/"

    type_a=["SNIa","SNIb","SNIc","SNII","SNIIn","SLSN-I","SLSN-II"]
    type_b=type_a

    allfiles=[]
    for sntype in sntypes:
        allfiles.extend(glob.glob(basedir+sntype+"/*.csv",recursive=True))
        #### DON'T CHANGE THIS
        #### Holdover from when the code output CSV rather than fixed-width ASCII
        #### Remaking the template files in the new format would take a week with
        ####   8 cores in parallel

    #type_b.sort()
    supscore = []; names=[' '] + type_b
    supscore = {}
    for itype in type_a:
        supscore[itype] = {}
        for jtype in type_b:
            if jtype == itype:
                continue
            supscore[itype][jtype] = {}
            supscore[itype][jtype]['scores'] = []
            supscore[itype][jtype]['phases'] = []
            supscore[itype][jtype]['truetypes'] = []
            supscore[itype][jtype]['names'] = []

    for file in allfiles:
        parts=file.split("/")
        stype=parts[-2]
        base=parts[-1]
        base=base.replace("sf_result_","")
        parts=base.split(".")
        snname=parts[0]
        ph=parts[1]
        truetype = None
        for itype in type_a+type_b:
            if itype in stype:
                truetype = itype
            if not(truetype):
                continue
        if ph=='max':
            phase=0.
        elif ph[0]=='m':
            phase=-1*float(ph[1:])
        elif ph[0]=='p':
            phase=1*float(ph[1:])
        tab=tasc.read(file)
        SN_matchs = tab['SN'].data

        delinds=[]
        for i,sn in enumerate(SN_matchs):
            if (snname in sn):
                delinds.append(i)
        SN_matchs=np.delete(SN_matchs, delinds)

        for itype in type_a:
            #if ("Ia" in itype) or ("Ic" in itype):
            #	wtyp1 = wheretype(itype, snfiles=SN_matchs, condense=True)
            #else:
            #	wtyp1 = wheretype(itype, snfiles=SN_matchs, type_trunc=data['temp_sn_tr'])
            wtyp1=wheretype(itype,SN_matchs,conIIn=True)
            ii1 = np.mean(wtyp1[:5])
            for jtype in type_b:
                if jtype == itype:
                    continue
                #if ("Ia" in jtype) or ("Ic" in jtype):
                #	wtyp2 = wheretype(jtype, snfiles=SN_matchs, condense=True)
                #else:
                #	wtyp2 = wheretype(jtype, snfiles=SN_matchs, type_trunc=data['temp_sn_tr'])
                wtyp2=wheretype(jtype,SN_matchs,conIIn=True)
                #dinds=[]
                #for i, ww in enumerate(wtyp2):
                #	if "IIb" in SN_matchs[ww]:
                #		dinds.append(i)
                #wtyp2=np.delete(wtyp2,dinds)
                ii2 = np.mean(wtyp2[:5])
                print(f"{itype}-{jtype}: {ii1-ii2:.2f}, phase={phase:.4f}, type={truetype}")
                supscore[itype][jtype]['scores'].append(ii1 - ii2)
                supscore[itype][jtype]['phases'].append(phase)
                supscore[itype][jtype]['truetypes'].append(truetype)
                supscore[itype][jtype]['names'].append(snname)
    with open(f"./tempbank/score_eval/scores_{basetype}.pickle", "wb") as f:
        pickle.dump(supscore, f)
    f.close()
    return


if __name__ == '__main__':

    sntypes=["SNIa-BL","SNIa-CL","SNIa-CN","SNIa-SS"]
    sntypes2=["SNIb","SNIb-pec"]
    sntypes3=["SNIc","SNIc-bl"]
    sntypes4=["SNII"]
    sntypes5=["SNIIn","SNIa-CSM","SLSN-IIn"]
    sntypes6=["SLSN-I"]
    sntypes7=["SLSN-II"]

    main(sntypes,"SNIa")
    main(sntypes2,"SNIb")
    main(sntypes3,"SNIc")
    main(sntypes4,"SNII")
    main(sntypes5,"SNIIn")
    main(sntypes6,"SLSN-I")
    main(sntypes7,"SLSN-II")
