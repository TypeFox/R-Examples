eval.penalty  <- function(basisobj, Lfdobj=int2Lfd(0),
                          rng=rangeval)
{
#  EVAL_PENALTY evaluates the inner products of a linear
#  differential operator L defined by LFDOBJ applied to a set of
#  basis functions defined by BASISOBJ.
#
#  LFDOBJ is a functional data object defining the order m
#  NONHOMOGENEOUS linear differential operator of the form
#  Lx(t) = w_0(t) x(t) + ... + w_{m-1}(t) D^{m-1}x(t) +
#          \exp[w_m(t)] D^m x(t).
#  This is a change from previous usage where LFDOBJ was assumed to
#  define a HOMOGONEOUS differential operator.  See function
#  $Lfd/Lfd() for details.
#
#  Arguments:
#  BASISOBJ ... Either a basis object or an fd object or
#               an fdPar object.  If an fdPar object,
#               and if no LFDOBJ is supplied, the LFDOBJ
#               in the fdPar object is used.
#  LFDOBJ   ... A linear differential operator object
#               applied to the functions that are evaluated.
#  RNG      ... A range over which the product is evaluated
#
#  Returns:
#  PENALTYMAT ... Symmetric matrix containing inner products.
#                 This matrix should be non-negative definite
#                 With NDERIV zero eigenvalues, where NDERIV
#                 is the highest order derivative in LFDOBJ.
#                 However, rounding error will likely cause
#                 NDERIV smallest eigenvalues to be nonzero,
#                 so be careful about calling CHOL or otherwise
#                 assuming the range is N - NDERIV.

#  last modified 21 December 2012

#  check BASISOBJ

if (inherits(basisobj, "fd")) basisobj <- basisobj$basis
	
if (inherits(basisobj, "fdPar")) {
    fdobj       <- basisobj$fd
    basisobj    <- fdobj$basis
}

if (!inherits(basisobj, "basisfd"))  stop(
	"Argument BASISOBJ is not a functional basis object.")

#  set up default values

rangeval <- basisobj$rangeval

#  deal with the case where LFDOBJ is an integer

Lfdobj <- int2Lfd(Lfdobj)

#  determine basis type

type <- basisobj$type

#  choose appropriate penalty matrix function

if (type=="bspline") penaltymat <- bsplinepen(basisobj, Lfdobj, rng)
else if(type=="const") {
        rangeval   <- getbasisrange(basisobj)
        penaltymat <- rangeval[2] - rangeval[1]
      }
else if(type=="expon")    penaltymat <- exponpen(basisobj,   Lfdobj)
else if(type=="fourier")  penaltymat <- fourierpen(basisobj, Lfdobj)
else if(type=="monom")    penaltymat <- monomialpen(basisobj,   Lfdobj)
else if(type=="polyg")    penaltymat <- polygpen(basisobj,   Lfdobj)
else if(type=="power")    penaltymat <- powerpen(basisobj,   Lfdobj)
else stop("Basis type not recognizable, can not find penalty matrix")

#  If drop indices are provided, drop rows and cols
#  associated with these indices

dropind <- basisobj$dropind
nbasis  <- basisobj$nbasis

if (length(dropind) > 0) {
    index <- 1:nbasis
    index <- index[-dropind]
    penaltymat <- penaltymat[index,index]
}

#  Make matrix symmetric since small rounding errors can
#  sometimes results in small asymmetries

penaltymat <- (penaltymat + t(penaltymat))/2

return(penaltymat)

}
