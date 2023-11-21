import numpy as np

def get_contset_placed(wav, cspacing, lines=None, zguess=0.):
    nwav = np.size(wav)
    placed = wav[0]
    for i in range(1,nwav):
        step = wav[i] - np.max(placed)
        if step >= cspacing[i]:
            if lines is not None:
                seps = np.abs(lines["LINE"]*(1+zguess) - wav[i]) - lines["MASKRAD"]
                nw = np.sum(seps <= 0)
                if nw >= 1:
                    continue
            placed = np.append(placed, wav[i])
    placed = np.delete(placed, [0])
    if placed[-1] == wav[-1]:
        placed = np.delete(placed, [-1])
    return placed
