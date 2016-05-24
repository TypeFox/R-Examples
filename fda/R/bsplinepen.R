bsplinepen <- function(basisobj, Lfdobj=2, rng=basisobj$rangeval,
                       returnMatrix=FALSE)
{

#  Computes the Bspline penalty matrix.
#  Arguments:
#  BASISOBJ    a basis.fd object of type "bspline"
#  LFDOBJ      a linear differential operator object.
#  RNG         a range over which the product is evaluate
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.

#  Returns the penalty matrix.

#  Last modified July 6, 2012 by Spencer Graves
#    for rng, rangeval, and params of class Date and POSIXct
#  Previously modified 9 May by Jim Ramsay

#  check BASISOBJ

  if (!(inherits(basisobj, "basisfd"))) stop(
    "First argument is not a basis.fd object.")

#  check basis type

  type <- basisobj$type
  if (type != "bspline") stop("basisobj not of type bspline")

#  check LFDOBJ

  Lfdobj <- int2Lfd(Lfdobj)

#  get basis information

  nbasis <- basisobj$nbasis
  Params <- basisobj$params

#  if there are no internal knots, use the monomial penalty

  if (length(Params) == 0) {
    basisobj      <- create.monomial.basis(rng, nbasis, 0:(nbasis-1))
    penaltymatrix <- monomialpen(basisobj, Lfdobj, rng)
    return(penaltymatrix)
  }

#  normal case:  PARAMS is not empty

# as.numeric
  Rangeval  <- basisobj$rangeval
  Rng <- rng
  op <- options(warn=-1)
  rng <- as.numeric(Rng)
  rangeval <- as.numeric(Rangeval)
  params <- as.numeric(Params)
  options(op)
# check
  nNA <- sum(is.na(rng))
  if(nNA>0)
    stop('as.numeric(rng) contains ', nNA,
         ' NA', c('', 's')[1+(nNA>1)],
         ';  class(rng) = ', class(Rng))
  nNAr <- sum(is.na(rangeval))
  if(nNAr>0)
    stop('as.numeric(rangeval) contains ', nNAr,
         ' NA', c('', 's')[1+(nNAr>1)],
         ';  class(rangeval) = ', class(Rangeval))
  nNAp <- sum(is.na(params))
  if(nNAp>0)
    stop('as.numeric(params) contains ', nNAp,
         ' NA', c('', 's')[1+(nNAp>1)],
         ';  class(params) = ', class(Params))

  breaks    <- c(rangeval[1],params,rangeval[2])  #  break points
  nbreaks   <- length(breaks)
  ninterval <- nbreaks - 1    #  number of intervals defined by breaks

#  check break values

  if (length(breaks) < 2)
    stop("The length of argument breaks is less than 2.")

#  Find the highest order derivative in LFD

  nderiv <- Lfdobj$nderiv

  norder <- nbasis - length(params)

#  check for order of derivative being equal or greater than
#  order of spline

  if (nderiv >= norder) {
    cat("\n")
    cat(paste(" Derivative of order", nderiv,
              "cannot be taken for B-spline of order", norder,"\n"))
    cat(" Probable cause is a value of the nbasis argument\n")
    cat(" in function create.basis.fd that is too small.\n")
    stop()
}

#  check for order of derivative being equal to order of spline
#  minus one, in which case following code won't work.

  if (nderiv > 0 && nderiv == norder - 1)
    stop(paste("Penalty matrix cannot be evaluated for derivative of order ",
               nderiv, " for B-splines of order ", norder))

#  special case where LFD is D^NDERIV and NDERIV = NORDER - 1

  bwtlist <- Lfdobj$bwtlist
  isintLfd <- TRUE
  if (nderiv > 0) {
    for (ideriv in 1:nderiv) {
      fdj <- bwtlist[[ideriv]]
      if (!is.null(fdj)) {
        if (any(fdj$coefs != 0)) {
            isintLfd <- FALSE
            break
        }
      }
    }
  }

  if (isintLfd && nderiv == norder - 1) {
    #  special case of nderiv = norder - 1
    halfseq    <- (breaks[2:nbreaks] + breaks[1:(nbreaks-1)])/2
    halfmat    <- bsplineS(halfseq, breaks, norder, nderiv, returnMatrix)
    brwidth    <- diff(breaks)
    penaltymat <- (t(halfmat) %*% diag(brwidth) %*% halfmat)
    return(penaltymat)
  }

#  look for knot multiplicities within the range

  intbreaks    <- c(rng[1], params, rng[2])
  index        <- intbreaks >= rng[1] & intbreaks <= rng[2]
  intbreaks    <- intbreaks[index]
  nintbreaks   <- length(intbreaks)
  uniquebreaks <- min(diff(intbreaks)) > 0

#  if LFD is D^NDERIV, and there are no break multiplicities,
#  use exact computation

  if (isintLfd && rng[1] == rangeval[1] &&
    uniquebreaks      && rng[2] == rangeval[2]) {

    #  Set up the knot sequence

    onesv <- matrix(1,1,norder)
    knots <- c(rangeval[1]*onesv, breaks[2:(nbreaks-1)], rangeval[2]*onesv)

    # Construct  the piecewise polynomial representation

    polyorder <- norder - nderiv
    ndegree   <- polyorder - 1
    prodorder <- 2*ndegree + 1   # order of product
    polycoef  <- array(0, c(ninterval, polyorder, norder))
    indxdown  <- seq(norder, (nderiv+1), -1)
    for (i in 1:nbasis) {
        #  compute polynomial representation of B(i,norder,t)(x)
        t       <- knots[i:(i+norder)]
        ppBlist <- ppBspline(t)
        Coeff   <- ppBlist[[1]]
        index   <- ppBlist[[2]]
        nrowcoef <- dim(Coeff)[1]
        # convert the index of the breaks in t to the index in the
        # variable "breaks"
        index <- index + i - norder
        CoeffD <- matrix(Coeff[,1:polyorder],nrowcoef,polyorder)
        if (nderiv > 0) {
            for (ideriv in 1:nderiv) {
                fac    <- indxdown - ideriv
                CoeffD <- outer(rep(1,nrowcoef),fac) * CoeffD
            }
        }
        # add the polynomial representation of B(i,norder,t)(x) to f
        if (i >= norder) k <- norder else k <- i
        if (i <= norder) m <- i      else m <- norder
        for (j in 1:nrowcoef) polycoef[i-k+j,,m-j+1] <- CoeffD[j,]
    }

    # Compute the scalar products

    prodmat <- matrix(0, nbasis, nbasis)
    convmat <- array(0,c(norder, norder, prodorder))
    for (k in 1:ninterval) {
        #  get the coefficients for the polynomials for this interval
        Coeff <- polycoef[k,,]
        #  compute the coefficients for the products
        for (i in 0:(ndegree-1)) {
            ind <- (0:i) + 1
            if (length(ind) == 1) {
                convmat[,,i+1        ] <-
                    outer(Coeff[ind,          ], Coeff[i-ind+2,      ])
                convmat[,,prodorder-i] <-
                    outer(Coeff[ndegree-ind+2,], Coeff[ndegree-i+ind,])
            } else {
                convmat[,,i+1        ] <-
                    crossprod(Coeff[ind,          ], Coeff[i-ind+2,      ])
                convmat[,,prodorder-i] <-
                    crossprod(Coeff[ndegree-ind+2,], Coeff[ndegree-i+ind,])
            }
        }
        ind <- (0:ndegree)+1
        convmat[,,ndegree+1] <-
                crossprod(Coeff[ind,          ], Coeff[ndegree-ind+2,])
        #  compute the coefficients of the integral
        delta    <- breaks[k+1] - breaks[k]
        power    <- delta
        prodmati <- matrix(0, norder, norder)
        for (i in 1:prodorder) {
            prodmati <- prodmati +
                power*convmat[,,prodorder-i+1]/i
            power <- power*delta
        }
        # add the integral to s
        index <- k:(k+norder-1)
        prodmat[index,index] <- prodmat[index,index] + prodmati
    }

     penaltymat <- prodmat

  } else {

    #  set iter

    iter <- 0

    #  LFDOBJ is not D^NDERIV, use approximate integration by calling
    #  function INPROD().

    if (uniquebreaks) {
        #  no knot multiplicities
        prodmat <- inprod(basisobj, basisobj, Lfdobj, Lfdobj, rng)
    } else {
        #  knot multiplicities  find their locations
        rngvec <- rng[1]
        for (i in 2:nbreaks) {
            if (breaks[i] == breaks[i-1]) rngvec <- c(rngvec, breaks[i])
        }
        rngvec <- unique(rngvec)
        nrng   <- length(rngvec)
        if (rngvec[nrng] < rng[2]) {
            rngvec <- c(rngvec,rng[2])
            nrng   <- nrng + 1
        }
        #  sum prodmat over intervals between knot multiplicities
        prodmat <- matrix(0,nbasis,nbasis)
        for (i in 2:nrng) {
            rngi <- c(rngvec[i-1] + 1e-10, rngvec[i] - 1e-10)
            prodmati <- inprod(basisobj, basisobj, Lfdobj, Lfdobj, rngi)
            prodmat <- prodmat + prodmati
        }
    }

    penaltymat <- prodmat
  }

    return( penaltymat )
}
