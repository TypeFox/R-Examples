OsplinePen=function(Boundary.knots, knots, ord=1)
{
  omega=if(ord==1){1
        }else if (ord==2){c(1,4,1)/3
        }else if (ord==3){c(14,64,24,64,14)/45
        }else if (ord==4){c(41,216,27,272,27,216,41)/140
        }else {
          k.seq=seq(0, 2*ord-2, by=1)
          omega=numeric(2*ord-1)
          for(k in 0:(2*ord-2)){
            f=function(t){
                  integrand=exp(sum(log(abs((t-k.seq)[-(k+1)])))-lfactorial(k)-lfactorial(2*ord-2-k))
                  integrand=integrand*prod(sign((t-k.seq)[-(k+1)]))
                  integrand[is.na(integrand)]=Inf
                  integrand}
            
            omega[k+1]=(-1)^k* integrate(Vectorize(f), 1e-10, 2*ord-2)$value
          }
          omega
        }
  Boundary.knots=sort(Boundary.knots)
  knots=sort(knots)
  K=length(knots)
  all.knots=c(rep(Boundary.knots[1L],2*ord), sort(knots), rep(Boundary.knots[2],2*ord))
  k.seq=seq(1, length(knots)+2*ord, by=1)
  h=diff(all.knots)/max(1,2*ord-2)
  w=xtilde=numeric((2*ord-1)*(K+4*ord-1))
  for(l in seq(1, K+4*ord-1, by=1))
    for(lp in seq(0, 2*ord-2, by=1)){
      xtilde[(2*ord-1)*(l-1)+lp+1]=all.knots[l]+lp*h[l]
      w[(2*ord-1)*(l-1)+lp+1]=h[l]*omega[lp+1]
    }
  Bder=splineDesign(all.knots, xtilde, ord=2*ord, derivs=rep(ord, length(xtilde)))
  Omega=crossprod(Bder, w*Bder)  
}

if(FALSE){
    library(splines)
    b.k=c(0,1)
    br=seq(.1,.9,by=.1)
    #    debug(OsplinePen)
    O1=OsplinePen(b.k, br, 1)
    O2=OsplinePen(b.k, br, 2)
    O3=OsplinePen(b.k, br, 3)
    O4=OsplinePen(b.k, br, 4)
    O5=OsplinePen(b.k, br, 5)
    O6=OsplinePen(b.k, br, 6)

    library(fda)
    des1=create.bspline.basis(c(0,1),norder=2, breaks=c(b.k[1], br, b.k[2]))
    P1=bsplinepen(des1, 1)
    P1=bsplinepen.fda(des1,1)
    max(abs(P1-O1))

    des2=create.bspline.basis(c(0,1),norder=4, breaks=c(b.k[1], br, b.k[2]))
    P2=bsplinepen(des2, 2)
    max(abs(P2-O2))

    des3=create.bspline.basis(c(0,1),norder=6, breaks=c(b.k[1], br, b.k[2]))
    P3=bsplinepen(des3, 3)
    max(abs(P3-O3))

    des4=create.bspline.basis(c(0,1),norder=8, breaks=c(b.k[1], br, b.k[2]))
    P4=bsplinepen(des4, 4, c(0,1))
    max(abs((P4-O4)/(P4+O4)*2),na.rm=T)

    des5=create.bspline.basis(c(0,1),norder=10, breaks=c(b.k[1], br, b.k[2]))
    P5=bsplinepen(des5, 5, c(0,1))
    max(abs((P5-O5)/(P5+O5)*2),na.rm=T)

    des6=create.bspline.basis(c(0,1),norder=12, breaks=c(b.k[1], br, b.k[2]))
    P6=bsplinepen(des6, 6, c(0,1))
    max(abs((P6-O6)/(P6+O6)*2),na.rm=T)


}

bsplinepen.fda <- function(basisobj, Lfdobj=2, rng=basisobj$rangeval)
{
    ### NOTE BY LQ: this function is the same as fda::bsplinepen, except that 
    ### the code to does not report error when nderiv == norder - 1 

    #  Computes the Bspline penalty matrix.
    #  Arguments:
    #  BASISOBJ    a basis.fd object of type "bspline"
    #  LFDOBJ      a linear differential operator object.
    #  RNG         a range over which the product is evaluate
    #  Returns the penalty matrix.

    #  Last modified 2008.08.23 by Spencer Graves
    #  Previously modified 26 October 2005

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
    params <- basisobj$params

    #  if there are no internal knots, use the monomial penalty

    if (length(params) == 0) {
        basisobj      <- create.monomial.basis(rng, nbasis, 0:(nbasis-1))
        penaltymatrix <- monomialpen(basisobj, Lfdobj, rng)
        return(penaltymatrix)
    }

    #  normal case:  PARAMS is not empty

    rangeval  <- basisobj$rangeval
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
        cat(paste("Derivative of order", nderiv,
                      "cannot be taken for B-spline of order", norder,"\n"))
        cat("Probable cause is a value of the nbasis argument\n")
        cat(" in function create.basis.fd that is too small.\n")
        stop()
    }

    #  check for order of derivative being equal to order of spline
    #  minus one, in which case following code won't work.

#    ### COMMENTED OUT BY LQ
#    if (nderiv > 0 && nderiv == norder - 1)
#        stop(paste("Penalty matrix cannot be evaluated for derivative of order ",
#                   nderiv, " for B-splines of order ", norder))
#    ### END OF COMMENT BY BY LQ

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
        halfmat    <- bsplineS(halfseq, breaks, norder, nderiv)
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
