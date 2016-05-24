# Function from package 'fda' (c) 2014

getbasismatrix <- function(evalarg, basisobj, nderiv=0, returnMatrix=FALSE) {
#  Computes the basis matrix evaluated at arguments in EVALARG associated
#    with basis.fd object BASISOBJ.  The basis matrix contains the values
#    at argument value vector EVALARG of applying the nonhomogeneous
#    linear differential operator LFD to the basis functions.  By default
#    LFD is 0, and the basis functions are simply evaluated at argument
#    values in EVALARG.
#
#  If LFD is a functional data object with m + 1 functions c_1, ... c_{m+1},
#    then it is assumed to define the order m HOMOGENEOUS linear
#    differential operator
#  Lx(t) = c_1(t) + c_2(t)x(t) + c_3(t)Dx(t) + ... +
#                             c_{m+1}D^{m-1}x(t) + D^m x(t).
#
#  If the basis type is either polygonal or constant, LFD is ignored.
#
#  Arguments:
#  EVALARG...Either a vector of values at which all functions are to evaluated,
#            or a matrix of values, with number of columns corresponding to
#            number of functions in argument FD.  If the number of evaluation
#            values varies from curve to curve, pad out unwanted positions in
#            each column with NA.  The number of rows is equal to the maximum
#            of number of evaluation points.
#  BASISOBJ...A basis object
#  NDERIV ... A nonnegative integer indicating a derivative to be evaluated.
#  RETURNMATRIX...If False, a matrix in sparse storage model can be returned
#             from a call to function BsplineS.  See this function for
#             enabling this option.

#
#  Note that the first two arguments may be interchanged.
#
#  Last modified 19 March 2014

##
##  Exchange the first two arguments if the first is an BASIS.FD object
##    and the second numeric
##
  if (is.numeric(basisobj) && inherits(evalarg, "basisfd")) {
    temp     <- basisobj
    basisobj <- evalarg
    evalarg  <- temp
  }
##
##  check EVALARG
##
#  if (!(is.numeric(evalarg)))  stop("Argument EVALARG is not numeric.")
  if(is.null(evalarg)) stop('evalarg required;  is NULL.')
  Evalarg <- evalarg
# turn off warnings in checking if argvals can be converted to numeric
  op <- options(warn=-1)
  evalarg <- as.numeric(Evalarg)
  options(op)
  nNA <- sum(is.na(evalarg))
  if(nNA>0)
    stop('as.numeric(evalarg) contains ', nNA,
         ' NA', c('', 's')[1+(nNA>1)],
         ';  class(evalarg) = ', class(Evalarg))

#  check basisobj

  if (!(inherits(basisobj, "basisfd")))
      stop("Second argument is not a basis object.")

#  search for stored basis matrix

  if (!(length(basisobj$basisvalues) == 0 || is.null(basisobj$basisvalues))) {
    #  one or more stored basis matrices found,
    #  check that requested derivative is available
    if (!is.vector(basisvalues)) stop("BASISVALUES is not a vector.")
    basisvalues <- basisobj$basisvalues
    nvalues     <- length(basisvalues)
    #  search for argvals match
    N  <- length(evalarg)
    OK <- FALSE
    for (ivalues in 1:nvalues) {
      basisvaluesi <- basisvalues[ivalues]
      if (!is.list(basisvaluesi)) stop("BASISVALUES does not contain lists.")
      argvals <- basisvaluesi[[1]]
      if (!length(basisvaluesi) < nderiv+2) {
          if (N == length(argvals)) {
              if (all(argvals == evalarg)) {
                  basismat <- basisvaluesi[[nderiv+2]]
                  OK <- TRUE
              }
          }
      }
    }
#   dimnames
    dimnames(basismat) <- list(NULL, basisobj$names)

    if (OK){
        if((!returnMatrix) && (length(dim(basismat)) == 2)){
            return(as.matrix(basismat))
        }
        return(basismat)
    }
  }

#  Extract information about the basis

  type     <- basisobj$type
  nbasis   <- basisobj$nbasis
  params   <- basisobj$params
  rangeval <- basisobj$rangeval
  dropind  <- basisobj$dropind

#  -----------------------------  B-spline basis  -------------------

  if (type == "bspline") {
      if (length(params) == 0) {
          breaks   <- c(rangeval[1], rangeval[2])
      } else {
   	    breaks   <- c(rangeval[1], params, rangeval[2])
      }
   	norder   <- nbasis - length(breaks) + 2
   	basismat <- bsplineS(evalarg, breaks, norder, nderiv, returnMatrix)

#  -----------------------------  Constant basis  --------------------

  } else if (type == "const") {
   	basismat  <- matrix(1,length(evalarg),1)

#  -------------------------------  Fourier basis  -------------------

  } else if (type == "fourier") {
   	period   <- params[1]
   	basismat <- fourier(evalarg, nbasis, period, nderiv)

#  -----------------------  Unrecognizable basis  --------------------

  } else {
   	stop("Basis type not recognizable")
  }
#  dimnames
  dimnames(basismat) <- list(NULL, basisobj$names)

#  remove columns for bases to be dropped

  if (length(dropind) > 0) basismat <- basismat[,-dropind]
  if (length(evalarg) == 1) {
    basismat = matrix(basismat,1,length(basismat))
  }
    

  if((!returnMatrix) && (length(dim(basismat)) == 2)){
    #  coerce basismat to be nonsparse
      return(as.matrix(basismat))
  } else {
    #  allow basismat to be sparse if it already is
      return(as.matrix(basismat))
  }

}
