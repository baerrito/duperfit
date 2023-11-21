"""
    MODULE: BSPLINE_FITTING
    *************************************
    Contains functions for quick B-Spline fitting & iterative fitting
"""

import numpy as np
from scipy.interpolate import splrep, BSpline

# ********************************************
# Default knot-setting function
# ********************************************

def default_knotfun(xdata, kspace=1.):
    """
    FUNCTION: DEFAULT_KNOTFUN
    ****************************************************************************
    Computes an array of knot-points for a given spacing and xdata. May not be
    optimal for some datasets.
    ****************************************************************************
    INPUT:
        xdata   :   Array along which to define knotpoints | real
        kspace  :   Spacing between knot points | real, default 1
    OUTPUT:
        knts    :   Array containing knotpoints, sans xmin and xmax, for use
                        in splrep | real
    ****************************************************************************
    Written by Michael Baer
    """
    xrange = xdata.max() - xdata.min()
    nkpts = int(xrange/float(kspace)) + 1
    tempkspace = xrange/float(nkpts-1)
    inds = np.array([i for i in range(0,nkpts)])
    knts = inds*tempkspace + xdata.min()
    knts = np.delete(knts, [0,-1])
    return knts

# ********************************************
# B-Spline fitting
# ********************************************

def bspline_fit(xdata, ydata, yerr=None, ivar=None, knots=None, knotfun=None,
                mask=None, **kwargs):
    """
    FUNCTION: BSPLINE_FIT
    ****************************************************************************
    Utilizes provided abscissa and ordinate points, either inverse variance or
    error in ordinate, either pre-defined knot points or a user-provided function
    for defining knot points, and any array masking to compute a weighted
    least-squares B-Spline fit.
    ****************************************************************************
    INPUT:
        xdata   :   Abscissa | real
        ydata   :   Ordinate | real
      At least one of the following is required:
        yerr    :   Error in ordinate | real, default None
        ivar    :   Inverse variance in ordinate | real, default None
      At least one of the following is required:
        knotfun :   Function for defining knot points | real function, default
                        None
        knots   :   Array containing pre-defined knot points without xdata[0]
                        and xdata[-1] | real, default None
    OPTIONAL INPUT:
        mask    :   Boolean array for any needed data masking | bool, default
                        None
        kwargs  :   Any necessary keyword arguments for knotfun
    OUTPUT:
        spl     :   B-Spline representation | real function
        outerr  :   Error code; =0, all is good
                                =1, no error or inverse variance provided
                                =2, no method for defining knot points provided
    ****************************************************************************
    Written by Michael Baer
    """
    if ivar is None:
        if yerr is None:
            print("BSPLINE_FIT: Error")
            print("Provide either error or inverse variance.")
            spl = 0.
            outerr=1
            return
        ivar = (yerr)**(-2.0)
    if mask is not None:
        xfit = xdata[mask]
        yfit = ydata[mask]
        ivar = ivar[mask]
    else:
        xfit = xdata
        yfit = ydata
    xsort = np.argsort(xdata)
    if knots is None:
        knots = knots
    if knotfun is not None:
        knots = knotfun(xdata[xsort], **kwargs)
    else:
        print("BSPLINE_FIT: Error")
        print("Give either a function for defining knot points")
        print("or an array of pre-defined knot points.")
        spl = 0.
        outerr = 2
        tck = 0.
        return spl, outerr, tck
    xsort = np.argsort(xfit)
    try:
        tck = splrep(xfit[xsort], yfit[xsort], w=np.sqrt(ivar[xsort]), k=3, task=-1, t=knots)
        spl = BSpline(*tck)
        outerr=0
    except:
        spl = 0.
        tck = 0.
        outerr = 3
    #tck = splrep(xfit[xsort], yfit[xsort], w=np.sqrt(ivar[xsort]), k=3, task=-1, t=knots)
    #spl = BSpline(*tck)
    #outerr=0
    return spl, outerr, tck

# ********************************************
# B-Spline iterative fitting with rejection
# ********************************************

def bspline_iterfit(xdata, ydata, yerr=None, ivar=None, knotfun=None, upper=5, lower=5,
                    maxiter=10, outmask=None, silence=False, **kwargs):
    """
    FUNCTION: BSPLINE_ITERFIT
    ****************************************************************************
    Utilizes provided abscissa and ordinate points, either inverse variance or
    error in ordinate, either pre-defined knot points or a user-provided
    function for defining knot points, and any array masking to compute a
    weighted least-squares B-Spline fit, iterating with rejection.
    ****************************************************************************
    INPUT:
        xdata   :   Abscissa | real
        ydata   :   Ordinate | real
      At least one of the following is required:
        knotfun :   Function for defining knot points | real function, default
                        None
        knots   :   Array containing pre-defined knot points without xdata[0]
                        and xdata[-1] | real, default None
    OPTIONAL INPUT:
        yerr    :   Error in ordinate | real, default standard deviation of
                        abscissa
        ivar    :   Inverse variance in ordinate | real, default None
        upper   :   Upper multiple of sigma for rejection | real, default 5
        lower   :   Lower multiple of sigma for rejection | real, default 5
        maxiter :   Maximum number of iterations | integer, default 10
        outmask :   Boolean array for any needed data masking | bool, default
                        True for all data points
        kwargs  :   Any necessary keyword arguments for knotfun
    OUTPUT:
        cspl    :   B-Spline representation after iteration is complete | real
                        function
        outmask :   Boolean array for any needed data masking post-rejection |
                        bool
    ****************************************************************************
    Written by Michael Baer
    """
    # Compute inverse variance if not provided
    if ivar is None:
        # If no y-error provided, default to ordinate variance
        if yerr is None:
            # Initialize the arrays to be the same size as ordinate
            ivar = np.ones(np.size(ydata))
            yerr = np.ones(np.size(ydata))
            # Make y-errors equal to ordinate's standard deviation
            yerr *= np.std(ydata)
            # Make inverse variance equal to inverse of ordinate variance
            var = np.var(ydata)
            ivar /= var
        else:
            ivar = (yerr)**(-2.0)
    else:
        if yerr is None:
            yerr = 1./np.sqrt(ivar)
    # Sort the data in ascending order
    xsort = np.argsort(xdata)
    xdata = xdata[xsort]
    ydata = ydata[xsort]
    ivar = ivar[xsort]
    # Default mask to True at all points
    if outmask is None:
        outmask = np.ones(np.size(ydata))==1
    # Compute working mask array
    maskwork = (outmask*(ivar>0.))[xsort]
    these = np.where(maskwork)
    nthese = np.sum(maskwork==1)
    iiter = 0
    qdone = False
    # Since we are masking points, we set a knot point defining function
    # If none is provided, use the built-in default
    if knotfun is None:
        knotfun = default_knotfun
    # Iterate the B-Spline fit
    while iiter <= maxiter and not qdone:
        if not silence:
            print(str(iiter) + '...', end=' ')
        # Fit the B-Spline
        cspl, outerr, tck = bspline_fit(xdata, ydata, ivar=ivar, knotfun=knotfun, mask=maskwork,
                                   **kwargs)
        if outerr == 3:
            iiter += 1
            continue
        elif outerr != 0:
            break
            return
        inmask = maskwork
        ymodel = cspl(xdata)
        res = ydata - ymodel
        badness = 0.0*outmask
        qbad = res<(-lower*yerr)
        badness += (((-res/(yerr + (yerr == 0))) > 0) * qbad)
        qbad = res>(upper*yerr)
        badness += (((res/(yerr + (yerr == 0))) > 0) * qbad)
        badness *= inmask
        newmask = (badness == 0)
        newmask = (newmask & inmask)
        # If the mask no longer changes, end the iteration early
        qdone = np.sum(newmask != maskwork) == 0
        # Iterate masks and iteration count
        maskwork = newmask
        iiter += 1
    outmask[xsort] = maskwork
    cspl, outerr, tck = bspline_fit(xdata, ydata, yerr, knotfun=knotfun, mask=outmask,
                               **kwargs)
    if outerr != 0:
        return 0., outerr, 0.
    if not silence:
        print("Done!")
    return cspl, outmask, tck
