# Function from package 'fda' (c) 2014

#  setClass statement for "Lfd" class

# setClass("Lfd",
# 	representation(call="call", nderiv="integer", bwtlist="list"))

#  Generator function for class Lfd

Lfd = function(nderiv=0, bwtlist=vector("list",0))
{

#  LFD creates a linear differential operator object of the form
#
#  Lx(t) = b_0(t) x(t) +  + b_{m-1}(t) D^{m-1}x(t) + D^m x(t)
#  or
#  Lx(t) = b_0(t) x(t) +  + b_{m-1}(t) D^{m-1}x(t) + D^m x(t)
#          \exp[b_m(t) D^m x(t).
#
#  Function x(t) is operated on by this operator L, and the operator
#  computes a linear combination of the function and its first m
#  derivatives.   The function x(t) must be scalar.
#  The operator L is called the HOMOGENEOUS because it does
#  not involve input or forcing functions.
#
#  The linear combination of derivatives is defined by the weight
#  or coefficient functions b_j(t), and these are assumed to vary
#  over t, although of course they may also be constant as a
#  special case.  Each of these must be scalar function.
#
#  The weight coefficient for D^m is special in that it must
#  be positive to properly identify the operator.  This is why
#  it is exponentiated.  In most situations, it will be 0,
#  implying a weight of one, and this is the default.
#
#  Some important functions also have the capability of allowing
#  the argument that is an LFD object be an integer. They convert
#  the integer internally to an LFD object by INT2LFD().  These are:
#     EVAL.FD()
#     EVAL.MON()
#     EVAL.POS()
#     EVAL.BASIS()
#     EVAL.PENALTY()
#
#  Arguments:
#
#  NDERIV  the order of the operator, that is,
#          the highest order of derivative.  This corresponds
#          to m in the above equation.
#  BWTLIST  A list object of length either NDERIV or NDERIV + 1.
#          Each list contains an FD object, an FDPAR object or
#          a numerical scalar constant.
#          If there are NDERIV functions, then the coefficient of D^m
#          is set to 1 otherwise, function NDERIV+1 contains a function
#          that is exponentiated to define the actual coefficient.
#          bwtlist may also be a vector of length NDERIV or NDERIVE + 1
#          containing constants.  If this is the case, the constants
#          are used as coefficients for a constant basis function.
#          The default is a row vector of NDERIV zeros.
#
#  Returns:
#
#  LFDOBJ  a functional data object

# last modified 2007 May 3 by Spencer Graves
#  Previously modified 9 October 2005

#  check nderiv

if (!is.numeric(nderiv))
    stop("Order of operator is not numeric.")
if (nderiv != round(nderiv))
    stop("Order of operator is not an integer.")
if (nderiv < 0)
    stop("Order of operator is negative.")

#  check that bwtlist is either a list or a fd object

if (!inherits(bwtlist, "list") && !inherits(bwtlist, "fd") &&
    !is.null(bwtlist) && !missing(bwtlist))
	stop("BWTLIST is neither a LIST or a FD object")

#  if bwtlist is missing or NULL, convert it to a constant basis FD object

if (is.null(bwtlist)) {
   bwtlist <- vector("list", nderiv)
	if (nderiv > 0) {
        conbasis <- create.constant.basis()
	    for (j in 1:nderiv) bwtlist[[j]]  <- fd(0, conbasis)
    }
}

#  if BWTLIST is a fd object, convert to a list object.

if (inherits(bwtlist, "fd")) bwtlist <- fd2list(bwtlist)

#  check size of bwtlist

nbwt <- length(bwtlist)

if (nbwt != nderiv & nbwt != nderiv + 1)
    stop("The size of bwtlist inconsistent with NDERIV.")

#  check individual list entries for class
#  and find a default range

if (nderiv > 0) {
    rangevec <- c(0,1)
    for (j in 1:nbwt) {
        bfdj <- bwtlist[[j]]
        if (inherits(bfdj, "fdPar")) {
			  bfdj <- bfdj$fd
			  bwtlist[[j]] <- bfdj
		 }
        if (!inherits(bfdj, "fd") && !inherits(bfdj, "integer"))
            stop(paste("An element of BWTLIST contains something other ",
               " than an fd object or an integer"))
        if (inherits(bfdj, "fd")) {
	        bbasis   <- bfdj$basis
	        rangevec <- bbasis$rangeval
        } else {
            if (length(bfdj) == 1) {
                bwtfd <- fd(bfdj, conbasis)
                bwtlist[[j]] <- bwtfd
            }
            else stop("An element of BWTLIST contains a more than one integer.")
        }
    }

    #  check that the ranges are compatible

    for (j in 1:nbwt) {
        bfdj    <- bwtlist[[j]]
        if (inherits(bfdj, "fdPar")) bfdj <- bfdj$fd
        bbasis <- bfdj$basis
        btype  <- bbasis$type
        #  constant basis can have any range
        if (!btype == "const") {
            brange = bbasis$rangeval
            if (any(rangevec != brange)) stop(
                "Ranges are not compatible.")
        }
    }
}

#  Save call

Lfd.call <- match.call()

#  set up the Lfd object

#  S4 definition
# Lfdobj <- new("Lfd", call=Lfd.call, nderiv=nderiv, bwtlist=bwtlist)

#  S3 definition

Lfdobj <- list(call=Lfd.call, nderiv=nderiv, bwtlist=bwtlist)
oldClass(Lfdobj) <- "Lfd"

Lfdobj

}

