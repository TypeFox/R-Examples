smooth.basisPar <- function(argvals, y, fdobj=NULL, Lfdobj=NULL,
      lambda=0, estimate=TRUE, penmat=NULL,
      wtvec=NULL, fdnames=NULL, covariates=NULL,
                         method="chol", dfscale=1 ) {
#  This function acts as a wrapper for those who don't want to take the
#  step of setting up a functional parameter object before invoking function
#  smooth.basis.  It simply does this setup for the user.   See the help
#  file for smooth.basis for further details.
#
#  However, smooth.basisPar also sets up a default basis in the event that
#  argument fdobj is either NULL
#  (order 4 b-spline basis with breaks = argvals) or
#  a positive integer (same as the NULL case, but with order = integer).

##
## 1.  check fdobj
##

    if (is.null(fdobj)) {
      #  if fdobj is NULL, create an order 4 bspline basis with breaks equal to
      #  argvals (see help file for create.bspline.basis(argvals) for further
      #  details)
      fdobj <- create.bspline.basis(argvals)
    } else {
      if (is.numeric(fdobj)) {
        # if fdobj is a positive integer, use this integer as the order
        # of a bspline basis with breaks equal to argvals
        if (length(fdobj)==1) {
          if (round(fdobj) != fdobj)
            stop("'fdobj' is numeric but not an integer")
          if (fdobj <= 0)
            stop("'fdobj' is not a positive integer")
          fdobj <- create.bspline.basis(argvals, norder=fdobj)
        } else {
          #  fdobj is neither NULL nor numeric, so use whatever it is as an
          #  argument for function fd, which in turn will do further checking
          fdobj <- fd(fdobj)
        }
      }
    }
##
## 2.  fdPar:  set up the functional parameter object from arguments
##
  fdP <- fdPar(fdobj, Lfdobj=Lfdobj, lambda=lambda,
               estimate=estimate, penmat=penmat)
##
## 3.  smooth.basis:  carry out smoothing by a call to smooth.basis and
##     return the smoothList object that this function returns
##
  smooth.basis(argvals, y, fdP, wtvec=wtvec, fdnames=fdnames,
               covariates=covariates, method="chol", dfscale=dfscale)
}
