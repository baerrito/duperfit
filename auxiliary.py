"""
    MODULE: AUXILIARY
    *************************************
    Contains auxiliary functions for the code
"""

import numpy as np

def load_lines_add(lines, name, restwav, isabs, isfit=False):
    # Check if absorption line
    if isabs:
        ew = -1.0
    else:
        ew = 1.0
    # Append line as structured array
    dtype = [('ION', '<U12'), ('FIT', bool), ('LINE', float), ('EW', float),
             ('TIEDTO', float), ('TIEDRATIO', float), ('ISABS', bool), ('MASKRAD', float)]
    line = np.array( (name, isfit, restwav, ew, 0.0, 0.0, isabs, float(restwav)/500.0),
                   dtype=dtype)
    lines.append(line)
    return

def load_lines():

    lines = []

    # Absorption lines
    load_lines_add(lines, 'Al II', 1670.80, True)
    load_lines_add(lines, 'Al III', 1854.70, True)
    load_lines_add(lines, 'Al III', 1862.80, True)
    load_lines_add(lines, 'Cr II', 2056.30, True)
    load_lines_add(lines, 'Cr II', 2066.20, True)
    load_lines_add(lines, 'Fe II', 2344.20, True)
    load_lines_add(lines, 'Fe II', 2373.70, True)
    load_lines_add(lines, 'Fe II', 2382.80, True)
    load_lines_add(lines, 'Fe II', 2586.70, True)
    load_lines_add(lines, 'Fe II', 2600.20, True)
    load_lines_add(lines, 'Mg II', 2795.528, True, isfit=True)
    load_lines_add(lines, 'Mg II', 2802.704, True, isfit=True)
    load_lines_add(lines, 'Mg I', 2852.127, True)
    load_lines_add(lines, 'Ca II', 3933.66, True)
    load_lines_add(lines, 'Ca II', 3968.47, True)
    load_lines_add(lines, 'Na I', 5891.94, True)

    # Emission lines
    load_lines_add(lines, 'O II', 3726.1, False, isfit=True)
    load_lines_add(lines, 'O II', 3728.8, False, isfit=True)
    load_lines_add(lines, 'H-alpha', 6562.79, False, isfit=True)
    load_lines_add(lines, 'H-beta', 4861.35, False, isfit=True)
    load_lines_add(lines, 'H-gamma', 4340.472, False, isfit=True)
    load_lines_add(lines, 'H-delta', 4101.734 , False)
    load_lines_add(lines, 'H-epsilon', 3970.075, False)
    load_lines_add(lines, 'H-8', 3889.064, False)
    load_lines_add(lines, 'H-9', 3835.397, False)
    load_lines_add(lines, 'H-10', 3797.909, False)
    load_lines_add(lines, 'H-11', 3770.633, False)
    load_lines_add(lines, 'H-12', 3750.151, False)
    load_lines_add(lines, 'H-13', 3734.369, False)
    load_lines_add(lines, 'H-14', 3721.946, False)
    load_lines_add(lines, 'H-15', 3711.978, False)
    load_lines_add(lines, 'H-16', 3703.859, False)
    load_lines_add(lines, 'O III', 4363.209 , False)
    load_lines_add(lines, 'O III', 4958.911 , False, isfit=True)
    load_lines_add(lines, 'O III', 5006.843 , False, isfit=True)
    load_lines_add(lines, 'Ne III', 3868.76, False)
    load_lines_add(lines, 'Ne III', 3967.47, False)
    load_lines_add(lines, 'N II', 6548.05 , False)
    load_lines_add(lines, 'N II', 6583.45, False)
    load_lines_add(lines, 'S II', 6717.0 , False)
    load_lines_add(lines, 'S II', 6731.3 , False)
    load_lines_add(lines, 'O II', 7320.1, False)
    load_lines_add(lines, 'O II', 7330.2, False)
    load_lines_add(lines, 'O I', 6300, False)
    load_lines_add(lines, 'O I', 6364, False)
    load_lines_add(lines, 'He I', 5015.678, False)
    load_lines_add(lines, 'He I', 5875.62, False)
    load_lines_add(lines, 'He I', 7065.19, False)
    load_lines_add(lines, 'He I', 6678.151, False)
    load_lines_add(lines, 'He I', 4471.4802, False)
    load_lines_add(lines, 'He I', 3888.648, False)
    load_lines_add(lines, 'He I', 4026.1914, False)
    load_lines_add(lines, 'He I', 10830.2, False)
    load_lines_add(lines, 'He II', 3203, False)
    load_lines_add(lines, 'He II', 4685.5, False)
    load_lines_add(lines, 'He II', 5411.52, False)
    load_lines_add(lines, 'He II', 6560.10, False)
    load_lines_add(lines, 'Pa-alpha', 18751.3, False)
    load_lines_add(lines, 'Pa-beta', 12818.072, False)
    load_lines_add(lines, 'Pa-gamma', 10938.17, False)
    load_lines_add(lines, 'Pa-delta', 10049.8, False)
    load_lines_add(lines, 'Pa-epsilon', 9546.2, False)
    load_lines_add(lines, 'Pa-8', 9229.7, False)
    load_lines_add(lines, 'Pa-9', 9015.3, False)
    load_lines_add(lines, 'Pa-10', 8862.89, False)
    load_lines_add(lines, 'Fe III', 4986, False)  ### check this
    load_lines_add(lines, 'S III', 6312, False)  ### check this
    load_lines_add(lines, 'Ar III', 7136, False)
    load_lines_add(lines, 'S III', 9069, False)
    load_lines_add(lines, 'S III', 9532, False)
    load_lines_add(lines, 'He I', 3188, False)

    return np.sort(lines,order=['LINE'])

def wheretype(sn_type,snfiles,conIIn=False,conflash=False):
	"""
	Assumes all spectrum files are within a type directory
	"""
	wtyp=[]
	for ii, sn in enumerate(snfiles):
		thedir=sn[:sn.rfind("/")]
		thetype=thedir[thedir.rfind("/")+1:]
		if sn_type=="SNIIn" and conIIn:
			if (thetype=="SNIIn") or (thetype=="SNIa-CSM") or (thetype=="SLSN-IIn"):
				wtyp.append(ii)
		elif conflash and sn_type=="flash":
			if "flash" in thetype:
				wtyp.append(ii)
		else:
			if ("SLSN" not in thetype) and ("-n" not in thetype) and \
                ("-CSM" not in thetype) and (not conflash and 'flash' not in thetype):
				vtype=thetype.split("-")[0]
			else:
				vtype=thetype
			if vtype==sn_type:
				wtyp.append(ii)
	return np.array(wtyp)
