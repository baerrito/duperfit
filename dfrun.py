"""
    SCRIPT: DFRUN
    *************************************
    A sample run of Duperfit, using the Type Ia SN 2023uta
    Uses the file `params.json` by default which the user may modify to fit
        their needs
    Else, the user may provide their own JSON-formatted file with input
        parameters following the same structure
"""

import os
import sys
import json

try:
    dfdir=os.environ['DFDIR']
except KeyError:
    raise Exception("$DFDIR not recognized")

sys.path.append(dfdir)

from dfclass import Duperfit

#####################################################################
# This script is still being developed; more features to come later #
#####################################################################

if __name__ == '__main__':

    ### checks for system argument
    if len(sys.argv)<2:
        fname="params.json"
    else:
        fname=sys.argv[1]
        if ".json" not in fname:
            fname=fname+'.json'
        if not os.path.isfile(fname):
            raise OSError(f'{fname} not found')

    with open(fname,"r") as f:
        params=json.load(f)

    DF = Duperfit(params)
    DF.fit(runRefit=True)
    DF.plot(0,saveFigure=True)
    DF.MIDscore(plotScores=True,evaluatedTypes=['SNIa','SNIb','SNIc','SNII','SNIIn','SLSN-I'])
