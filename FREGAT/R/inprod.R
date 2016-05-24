# Function from package 'fda' (c) 2014

inprod <- function(fdobj1, fdobj2=NULL, Lfdobj1=int2Lfd(0), Lfdobj2=int2Lfd(0),
                   rng = range1, wtfd = 0, returnMatrix=FALSE)
{

#  computes matrix of inner products of functions by numerical
#    integration using Romberg integration

#  Arguments:
#  FDOBJ1 and FDOBJ2    These may be either functional data or basis
#               function objects.  In the latter case, a functional
#               data object is created from a basis function object
#               by using the identity matrix as the coefficient matrix.
#               Both functional data objects must be univariate.
#               If inner products for multivariate objects are needed,
#               use a loop and call inprod(FDOBJ1[i],FDOBJ2[i]).
#     If FDOBJ2 is not provided or is NULL, it defaults to a function
#     having a constant basis and coefficient 1 for all replications.
#     This permits the evaluation of simple integrals of functional data
#     objects.
#  LFDOBJ1 and LFDOBJ2  order of derivatives for inner product for
#               FDOBJ1 and FDOBJ2, respectively, or functional data
#               objects defining linear differential operators
#  RNG    Limits of integration
#  WTFD   A functional data object defining a weight
#  JMAX   maximum number of allowable iterations
#  EPS    convergence criterion for relative stop
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.

#  Return:
#  A matrix of NREP1 by NREP2 of inner products for each possible pair
#  of functions.

#  Last modified 24 December 2012 by Jim Ramsay

#  Check FDOBJ1 and get no. replications and basis object

result1   <- fdchk(fdobj1)
nrep1     <- result1[[1]]
fdobj1    <- result1[[2]]
coef1     <- fdobj1$coefs
basisobj1 <- fdobj1$basis
type1     <- basisobj1$type
range1    <- basisobj1$rangeval

#  Default FDOBJ2 to a constant function, using a basis that matches
#  that of FDOBJ1 if possible.

if (is.null(fdobj2)) {
    tempfd    <- fdobj1
    tempbasis <- tempfd$basis
    temptype  <- tempbasis$type
    temprng   <- tempbasis$rangeval
    if (temptype == "bspline") {
        basis2 <- create.bspline.basis(temprng, 1, 1)
    } else {
        if (temptype == "fourier") basis2 <- create.fourier.basis(temprng, 1)
        else                       basis2 <- create.constant.basis(temprng)
    }
    fdobj2 <- fd(1,basis2)
}

#  Check FDOBJ2 and get no. replications and basis object

result2   <- fdchk(fdobj2)
nrep2     <- result2[[1]]
fdobj2    <- result2[[2]]
coef2     <- fdobj2$coefs
basisobj2 <- fdobj2$basis
type2     <- basisobj2$type
range2    <- basisobj2$rangeval

# check ranges

if (rng[1] < range1[1] || rng[2] > range1[2]) stop(
	 "Limits of integration are inadmissible.")

#  Call B-spline version if
#  [1] both functional data objects are univariate
#  [2] both bases are B-splines
#  (3) the two bases are identical
#  (4) both differential operators are integers
#  (5) there is no weight function
#  (6) RNG is equal to the range of the two bases.

if (inherits(fdobj1,"fd")            && 
    inherits(fdobj2,"fd")            &&
    type1 == "bspline"               && 
    type2 == "bspline"               &&
    is.eqbasis(basisobj1, basisobj2) &&
    is.integer(Lfdobj1)              && 
    is.integer(Lfdobj2)              &&
    length(basisobj1$dropind) == 0   &&
    length(basisobj1$dropind) == 0   &&
    wtfd == 0                        && all(rng == range1)) {

    inprodmat <- inprod.bspline(fdobj1, fdobj2,
                     Lfdobj1$nderiv, Lfdobj2$nderiv)
    return(inprodmat)
}

#  check LFDOBJ1 and LFDOBJ2

Lfdobj1 <- int2Lfd(Lfdobj1)
Lfdobj2 <- int2Lfd(Lfdobj2)

#  Else proceed with the use of the Romberg integration.

#  ------------------------------------------------------------
#  Now determine the number of subintervals within which the
#  numerical integration takes.  This is important if either
#  basis is a B-spline basis and has multiple knots at a
#  break point.
#  ------------------------------------------------------------

#  set iter

iter <- 0

# The default case, no multiplicities.

rngvec <- rng

#  check for any knot multiplicities in either argument

knotmult <- numeric(0)
if (type1 == "bspline") knotmult <- knotmultchk(basisobj1, knotmult)
if (type2 == "bspline") knotmult <- knotmultchk(basisobj2, knotmult)

#  Modify RNGVEC defining subinvervals if there are any
#  knot multiplicities.

if (length(knotmult) > 0) {
    knotmult <- sort(unique(knotmult))
    knotmult <- knotmult[knotmult > rng[1] && knotmult < rng[2]]
    rngvec   <- c(rng[1], knotmult, rng[2])
}

#  check for either coefficient array being zero
if ((all(c(coef1) == 0) || all(c(coef2) == 0)))
	return(matrix(0,nrep1,nrep2))

#  -----------------------------------------------------------------
#                   loop through sub-intervals
#  -----------------------------------------------------------------

#  Set constants controlling convergence tests

JMAX <- 15
JMIN <-  5
EPS  <- 1e-4

inprodmat <- matrix(0,nrep1,nrep2)

nrng <- length(rngvec)
for (irng  in  2:nrng) {
    rngi <- c(rngvec[irng-1],rngvec[irng])
    #  change range so as to avoid being exactly on
    #  multiple knot values
    if (irng > 2   ) rngi[1] <- rngi[1] + 1e-10
    if (irng < nrng) rngi[2] <- rngi[2] - 1e-10

    #  set up first iteration

    iter  <- 1
    width <- rngi[2] - rngi[1]
    JMAXP <- JMAX + 1
    h <- rep(1,JMAXP)
    h[2] <- 0.25
    s <- array(0,c(JMAXP,nrep1,nrep2))
    sdim <- length(dim(s))
    #  the first iteration uses just the endpoints
    fx1 <- eval.fd(rngi, fdobj1, Lfdobj1, returnMatrix)
    fx2 <- eval.fd(rngi, fdobj2, Lfdobj2, returnMatrix)
    #  multiply by values of weight function if necessary
    if (!is.numeric(wtfd)) {
        wtd <- eval.fd(rngi, wtfd, 0, returnMatrix)
        fx2 <- matrix(wtd,dim(wtd)[1],dim(fx2)[2]) * fx2
    }
    s[1,,] <- width*as.numeric(crossprod(fx1,fx2))/2
    tnm  <- 0.5
    iter <- 1

    #  now iterate to convergence

    for (iter in 2:JMAX) {
        tnm <- tnm*2
        if (iter == 2) {
            x <- mean(rngi)
        } else {
            del <- width/tnm
            x   <- seq(rngi[1]+del/2, rngi[2]-del/2, del)
        }
        fx1 <- eval.fd(x, fdobj1, Lfdobj1, returnMatrix)
        fx2 <- eval.fd(x, fdobj2, Lfdobj2, returnMatrix)
        if (!is.numeric(wtfd)) {
            wtd <- eval.fd(wtfd, x, 0, returnMatrix)
            fx2 <- matrix(wtd,dim(wtd)[1],dim(fx2)[2]) * fx2
        }
        chs <- width*as.numeric(crossprod(fx1,fx2))/tnm
        s[iter,,] <- (s[iter-1,,] + chs)/2
        if (iter >= 5) {
            ind <- (iter-4):iter
            ya <- s[ind,,]
            ya <- array(ya,c(5,nrep1,nrep2))
            xa <- h[ind]
            absxa <- abs(xa)
            absxamin <- min(absxa)
            ns <- min((1:length(absxa))[absxa == absxamin])
            cs <- ya
            ds <- ya
            y  <- ya[ns,,]
            ns <- ns - 1
            for (m in 1:4) {
                for (i in 1:(5-m)) {
                    ho      <- xa[i]
                    hp      <- xa[i+m]
                    w       <- (cs[i+1,,] - ds[i,,])/(ho - hp)
                    ds[i,,] <- hp*w
                    cs[i,,] <- ho*w
                }
                if (2*ns < 5-m) {
                    dy <- cs[ns+1,,]
                } else {
                    dy <- ds[ns,,]
                    ns <- ns - 1
                }
                y <- y + dy
            }
            ss     <- y
            errval <- max(abs(dy))
            ssqval <- max(abs(ss))
            if (all(ssqval > 0)) {
                crit <- errval/ssqval
            } else {
                crit <- errval
            }
            if (crit < EPS && iter >= JMIN) break
        }
        s[iter+1,,] <- s[iter,,]
        h[iter+1]   <- 0.25*h[iter]
        if (iter == JMAX) warning("Failure to converge.")
    }

    inprodmat <- inprodmat + ss

}

if((!returnMatrix) && (length(dim(inprodmat)) == 2)) {
    #  coerce inprodmat to be nonsparse
    return(as.matrix(inprodmat))
} else {
    #  allow inprodmat to be sparse if it already is
    return(inprodmat)
}

}

#  -------------------------------------------------------------------------------

fdchk <- function(fdobj) {

    #  check the class of FDOBJ and extract coefficient matrix

    if (inherits(fdobj, "fd")) coef  <- fdobj$coefs
    else
	  if (inherits(fdobj, "basisfd")) {
    	    coef  <- diag(rep(1,fdobj$nbasis - length(fdobj$dropind)))
    	    fdobj <- fd(coef, fdobj)
	  }
    else stop("FDOBJ is not an FD object.")

    #  extract the number of replications and basis object

    coefd <- dim(as.matrix(coef))
    if (length(coefd) > 2) stop("Functional data object must be univariate")
    nrep     <- coefd[2]
    basisobj <- fdobj$basis

    return(list(nrep, fdobj))

}

#  -------------------------------------------------------------------------------

knotmultchk <- function(basisobj, knotmult) {
    type <- basisobj$type
    if (type == "bspline") {
        # Look for knot multiplicities in first basis
        params  <- basisobj$params
        nparams <- length(params)
        if (nparams > 1) {
            for (i in 2:nparams) {
                if (params[i] == params[i-1]) {
                    knotmult <- c(knotmult, params[i])
                }
            }
        }
    }
    return(knotmult)
}


