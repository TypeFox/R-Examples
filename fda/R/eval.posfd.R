predict.posfd <- function(object, newdata=NULL, Lfdobj=0, 
                          returnMatrix=FALSE, ...){
#  Last modified  7 May 2012 by Jim Ramsay
  if (is.null(newdata)) newdata <- object$argvals
  evalPos <- eval.posfd(newdata, object$Wfdobj, Lfdobj, returnMatrix=FALSE)
#
  evalPos
}

fitted.posfd <- function(object, ...){
  predict(object)
}

residuals.posfd <- function(object, ...){
  pred <- predict(object)
  object$y-pred
}

eval.posfd <- function(evalarg, Wfdobj, Lfdobj=int2Lfd(0), returnMatrix=FALSE)
{
#  Evaluates a value or a derivative of a positive functional
#  data object.
#  A positive functional data object h  is = the form
#           h(x) = (exp Wfdobj)(x)
#  Note that the linear differential operator object LFDOBJ
#  MUST be an integer in the range 0 to 1.
#  Note that the first two arguments may be interchanged.
#
#  Arguments:
#  EVALARG ... A vector of values at which all functions are to
#              evaluated.
#  WFDOBJ  ... Functional data object.  It must define a single
#              functional data observation.
#  LFDOBJ  ... A linear differential operator object
#              applied to the functions that are evaluated.
#              Default is INT2LFD(0).
#  RETURNMATRIX ... If FALSE, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.
#
#  Returns:  An array of function values corresponding to the
#              argument values in EVALARG

#  Exchange the first two arguments if the first is an FD object
#    and the second numeric

if (is.numeric(Wfdobj) & inherits(class(evalarg), "fd")) {
    temp    <- Wfdobj
    Wfdobj  <- evalarg
    evalarg <- temp
}

#  Check the arguments

if (!is.numeric(evalarg)) stop(
	"Argument EVALARG is not numeric.")

evalarg <- as.vector(evalarg)

#  check WFDOBJ

if (!inherits(Wfdobj, "fd")) stop(
    "Argument WFDOBJ is not a functional data object.")

#  check LFDOBJ

Lfdobj = int2Lfd(Lfdobj)
if (!inherits(Lfdobj, "Lfd")) stop(
    "LFDOBJ is not linear differential operator object.")

nderiv = Lfdobj$nderiv

#  Extract information about the basis

basisobj <- Wfdobj$basis
nbasis   <- basisobj$nbasis
rangeval <- basisobj$rangeval
onerow   <- rep(1,nbasis)

#  Set up coefficient array for FD

coef  <- Wfdobj$coefs

#  Evaluate function values

index <- evalarg < rangeval[1]-1e-10
if (length(evalarg[index]) > 0) evalarg <- evalarg[!index]
index <- evalarg > rangeval[2]+1e-10
if (length(evalarg[index]) > 0) evalarg <- evalarg[!index]

basismat <- getbasismatrix(evalarg, basisobj, 0, returnMatrix)
fdvec    <- exp(basismat %*% coef)

#  If a differential operator has been defined in LFDOBJ, compute
#  the derivative values

if (nderiv > 0) {
     Lbasismat <- eval.basis(evalarg, basisobj, Lfdobj, returnMatrix)
     evalarray <- fdvec*(Lbasismat %*% coef)
} else evalarray <- fdvec

return(evalarray)

}

