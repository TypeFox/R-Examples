#  setClass for "fdPar"

# setClass("fdPar", representation(fd       = "fd",
#                                  Lfd      = "Lfd",
#                                  lambda   = "numeric",
#                                  estimate = "logical",
#                                  penmat   = "matrix",))

#  Generator function of class fdPar

fdPar <- function(fdobj=NULL, Lfdobj=NULL, lambda=0, estimate=TRUE,
                  penmat=NULL){

# Sets up a functional parameter object
#  Arguments:
#  FDOBJ    ... A functional data object.
#               The basis for this object is used to define
#               the functional parameter, or functional
#               parameters of FDOBJ has replications.
#               When an initial value is required for iterative
#               estimation of a functional parameter, the coefficients
#               will give the initial values for the iteration.
#  LFDOBJ   ... A linear differential operator value or a derivative
#               value for penalizing the roughness of the object.
#               By default, this is 0.
#  LAMBDA   ... The penalty parameter controlling the smoothness of
#               the estimated parameter.  By default this is 0.
#  ESTIMATE ... If nonzero, the parameter is estimated; if zero, the
#               parameter is held fixed at this value.
#               By default, this is 1.
#  PENMAT   ... The penalty matrix.
#               In repeated calls to SMOOTH_BASIS, if this is
#               saved, then the penalty does not need evaluating
#               repeatedly.  Don't use, though, if LFDOBJ or LAMBDA
#               are changed in the calculation.
#
#  An alternative argument list:
#  The first argument can also be a basis object.  In this case, an
#  FD object is set up with an empty coefficient matrix.
#  For many purposes, the coefficient array is either not needed, or
#  supplied later.
#
#  Return:
#  FDPAROBJ ... A functional parameter object

#  Last modified 22 December 2012 by Jim Ramsay

#  ----------------------------------------------------------------------
#                            Default fdPar objects
#  ----------------------------------------------------------------------

  if(!inherits(fdobj, 'fd')) {
    if (is.null(fdobj)) {
    #  fdPar called without arguments
      fdobj = fd()
    }  else {
      if (inherits(fdobj, "basisfd")) {
        #  if the first argument is a basis object, convert it to
        #  a default FD object with an empty coefficient matrix.
        nbasis  <- fdobj$nbasis
        dropind <- fdobj$dropind
        nbasis  <- nbasis - length(dropind)
        coefs   <- matrix(0,nbasis,nbasis)
        fdnames <- list('time', 'reps 1', 'values')
        if(!is.null(fdobj$names)){
          basisnames <- {
            if(length(dropind)>0)
              fdobj$names[-dropind]
            else
              fdobj$names
          }
          dimnames(coefs) <- list(basisnames, NULL)
          fdnames[[1]] <- basisnames
        }
        fdobj <- fd(coefs, fdobj, fdnames)
      }
      else if(is.numeric(fdobj))fdobj <- fd(fdobj)

      else stop("First argument is neither a functional data object ",
                "nor a basis object.")
    }
  }

#  ----------------------------------------------------------------------
#                            Check parameters
#  ----------------------------------------------------------------------

#  check Lfdobj

  {
    if (is.null(Lfdobj)) {
      if(fdobj$basis$type=='fourier'){
        rng <- fdobj$basis$rangeval
        Lfdobj <- vec2Lfd(c(0,(2*pi/diff(rng))^2,0), rng)
#        warning("Provding default Lfdobj = harmonic acceleration ",
#                "operator on c(", rng[1], ', ', rng[2],
#                ') = vec2Lfd(c(0,(2*pi/diff(rng))^2,0), rng);',
#                '  [default prior to fda 2.1.0:  int2Lfd(0)].')
      } else {
        norder <- {
          if (fdobj$basis$type=='bspline') norder.bspline(fdobj$basis)
          else 2
        }
        Lfdobj <- int2Lfd(max(0, norder-2))
      }
    }
    else
      Lfdobj <- int2Lfd(Lfdobj)
  }

  if (!inherits(Lfdobj, "Lfd"))
    stop("'Lfdobj' is not a linear differential operator object.")

#  check lambda

if (!is.numeric(lambda)) stop("Class of LAMBDA is not numeric.")
if (lambda < 0) stop("LAMBDA is negative.")

#  check estimate

if (!is.logical(estimate)) stop("Class of ESTIMATE is not logical.")

#  check penmat

if (!is.null(penmat)) {
    if (!is.numeric(penmat)) stop("PENMAT is not numeric.")
#    penmatsize <- size(penmat)
    penmatsize <- dim(penmat)
    if (any(penmatsize != nbasis)) stop("Dimensions of PENMAT are not correct.")
}

#  ----------------------------------------------------------------------
#                    set up the fdPar object
#  ----------------------------------------------------------------------

#  S4 definition
# fdParobj <- new("fdPar", fd=fdobj, Lfd=Lfdobj, lambda=lambda, estimate=estimate,
#                  penmat=penmat)

#  S3 definition

fdParobj <- list(fd=fdobj, Lfd=Lfdobj, lambda=lambda, estimate=estimate,
                 penmat=penmat)

oldClass(fdParobj) <- "fdPar"

fdParobj

}

#  ----------------------------------------------------------------------

#  "print" method for "fdPar"

print.fdPar <- function(x, ...)
{
  object <- x
  cat("Functional parameter object:\n\n")
      print("Functional data object:")
  print.fd(object$fd)
      print("Linear differential operator object:")
  print.Lfd(object$Lfd)
  cat(paste("\nSmoothing parameter =",object$lambda,"\n"))
  cat(paste("\nEstimation status =",object$estimate,"\n"))
      if (!is.null(object$penmat)) {
          print("Penalty matrix:")
          print(object$penmat)
      }
}

#  ----------------------------------------------------------------------

#  "summary" method for "fdPar"

summary.fdPar <- function(object, ...)
{
  cat("Functional parameter object:\n\n")
      print("Functional data object:")
  summary.fd(object$fd)
      print("Linear differential operator object:")
  summary.Lfd(object$Lfd)
  cat(paste("\nSmoothing parameter =",object$lambda,"\n"))
  cat(paste("\nEstimation status =",object$estimate,"\n"))
      if (!is.null(object$penmat))
          print(paste("Penalty matrix dimensions:",dim(object$penmat)))
}


