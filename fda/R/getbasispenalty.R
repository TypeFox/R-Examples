getbasispenalty <- function(basisobj, Lfdobj=NULL)
{
#  Computes the penaltymat matrix  associated with basis object basisobj.
#    This is defined in terms of a linear differential operator LFDOBJ.
#    The default for LFDOBJ depends on the nature of the basis.

#  Last modified 19 March 2014

#  check BASISOBJ

if (!(inherits(basisobj, "basisfd"))) stop(
    "First argument is not a basis object.")

type   <- basisobj$type
nbasis <- basisobj$nbasis

if        (type == "fourier") {
    if (is.null(Lfdobj)) Lfdobj <- 2
    penaltymat <- fourierpen(basisobj, Lfdobj)
} else if (type == "bspline") {
    norder <- basisobj$nbasis - length( basisobj$params )
    if (is.null(Lfdobj)) Lfdobj <- 2
    penaltymat <- bsplinepen(basisobj, Lfdobj)
} else if (type == "expon")   {
    if (is.null(Lfdobj)) Lfdobj <- 2
    penaltymat <- exponpen(basisobj, Lfdobj)
} else if (type == "polyg" | type == "polygonal")   {
    if (is.null(Lfdobj)) Lfdobj <- 1
    penaltymat <- polygpen(basisobj, Lfdobj)
} else if (type == "power")   {
    if (is.null(Lfdobj)) Lfdobj <- 2
    penaltymat <- powerpen(basisobj, Lfdobj)
} else if (type == "const")   {
    if (is.null(Lfdobj)) Lfdobj <- 0
    if (Lfdobj == 0) {
      penaltymat <- basisobj$rangeval[2] - basisobj$rangeval[1]
    } else {
      penaltymat <- 0
    }
} else {
    stop("Basis type not recognizable")
}

return(penaltymat)
}
