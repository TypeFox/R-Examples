"lmoms.bootbarvar" <-
function(x, nmom=6, covarinverse=TRUE, verbose=TRUE,
            force.exact=FALSE, nohatSIGMA=FALSE, nsim=500, bign=40, ...) {
   # Begin with the definition of several helper functions
   # that follow the nomenclature of Hutson and Ernst (2000)
   # [ Journal Royal Statistical Society B, 62(1), pp. 89-94 ]
   # note that in these functions that the incomplete beta function
   # is computed by the CDF of the beta distribution multiplied by
   # a complete beta function. So the notation or use in this algorithm
   # does not mimic the nomenclature of Hutson and Ernst (2000) in
   # this regard.
   "WJR" <-
   function(j,r,n) { #  equation 2.3
      tmp <- pbeta( j/n, r, n-r+1 ) - pbeta( (j-1)/n, r, n-r+1 )
      tmp <- tmp * r * choose(n, r) * beta(r,n-r+1)
      return(tmp)
   }

   "WIJRS" <-
   function(i,j,r,s,n) { #  unnumbered equation page 93 (near bottom)
      nCrs <- exp(lgamma(n+1) - lgamma(r-1+1) - lgamma(s-r-1+1) - lgamma(n-s+1))
      afunc <- function(k) {
         v <- choose(s-r-1, k)*(-1)^(s-r-1-k)/(s-k-1)
         v <- v*( (i/n)^(s-k-1) - ((i-1)/n)^(s-k-1) )
         b <- pbeta( j/n, k+1, n-s+1 ) - pbeta( (j-1)/n, k+1, n-s+1 )
         b <- b*beta(k+1, n-s+1)
         return(v*b)
      }
      wijrs <- sapply(0:(s-r-1), afunc)
      wijrs <- nCrs*sum(wijrs)
      return(wijrs)
   }

   "VJRS" <-
   function(j,r,s,n) { #  unnumbered equation page 93 (bottom)
      nCrs <- exp(lgamma(n+1) - lgamma(r-1+1) - lgamma(s-r-1+1) - lgamma(n-s+1))
      afunc <- function(k) {
         v  <- choose(s-r-1, k)*(-1)^(s-r-1-k)/(s-k-1)
         b1 <- pbeta( j/n, s, n-s+1 ) - pbeta( (j-1)/n, s, n-s+1 )
         b1 <- b1*beta(s, n-s+1)
         b2 <- pbeta( j/n, k+1, n-s+1 ) - pbeta( (j-1)/n, k+1, n-s+1 )
         b2 <- b2*beta(k+1, n-s+1)
         v  <- v * (b1 - ((j-1)/n)^(s-k-1) * b2)
         return(v)
      }
      wijrs <- sapply(0:(s-r-1), afunc)
      wijrs <- nCrs*sum(wijrs)
      return(wijrs)
   }

   "hatMUk" <-
   function(X, k=1) { #  equation 2.2
      # Exact bootstrap estimates of order statistic expectations if k=1
      n <- length(X)
      uhatrk <- sapply(1:n, function(r) {
                 tmp <- sapply(1:n, function(j) { WJR(j,r,n)*X[j]^k })
                         return (sum(tmp)) })
      return(uhatrk)
   }

   "hatSIGMAsq" <-
   function(X) { #  equation 3.1
      n <- length(X)
      hatmu <- hatMUk(X)
      sigmasq <- sapply(1:n, function(r) { tmp <- sapply(1:n, function(j) {
                        WJR(j,r,n)*(X[j] - hatmu[r])^2 })
                        return(sum(tmp)) })
      return(sigmasq)
   }

   "hatSIGMA" <-
   function(X) { #  equation 3.2
      n <- length(X)
      muvec1   <- hatMUk(X)
      varcovar <- matrix(nrow=n, ncol=n)
      diag(varcovar) <- hatSIGMAsq(X)
      if(verbose) message("Generating the hatSIGMA matrix (all else is reasonably fast)")
      for(r in 1:n) {
         if(verbose) message("r=",r,", s:", appendLF = FALSE)
         for(s in 1:n) {
            if(r >= s) next
            if(verbose) message(" ",s, appendLF = FALSE)
            tmpA <- sapply(2:n, function(j) {
               tmp <- sapply(1:(j-1),
                         function(i) { WIJRS(i,j,r,s,n)   *
                                       (X[i] - muvec1[r]) *
                                       (X[j] - muvec1[s]) } )
                                         return(sum(tmp)) } )
            tmpA <- sum(tmpA)
            tmpB <- sapply(1:n, function(j) { VJRS(j,r,s,n) *
                                              (X[j] - muvec1[r]) *
                                              (X[j] - muvec1[s]) } )
            tmpB <- sum(tmpB)
            z <- tmpA + tmpB
            if(is.na(z)) z <- 0
            varcovar[r,s] <- z
         }
         if(verbose) message("")
      }
      if(verbose) message("")
      varcovar[lower.tri(varcovar)] <- t(varcovar)[lower.tri(varcovar)]
      return(varcovar)
   }
   #  functional definitions are complete

    # Begin real logic
    x <- sort(x) #  generate the order statistics
    n <- length(x)
    if(nmom > n) {
       stop("More L-moments requested by parameter 'nmom' than data points available in 'x'")
    }
    if(nmom <= 2) {
       stop("Too few L-moments requested in nmom (r >= 3)")
    }
    if(length(unique(x)) == 1) stop("all values are equal--Lmoments can not be computed")
    lmr <- lmoms(x, nmom=nmom) # regular sample L-moments

    if(! force.exact & n > bign) {
       warning("Sample size is >bign, testing with highly non-normal data suggests this as an ",
               "upper limit for numerical stability, will proceed with simulations (nsim) and ",
               "estimated variances only, unless force.exact=TRUE, the return L-moments and ratios ",
               "simply come from lmoms(), a hint of numerical problems with the exact are negative variances")
       LA <- LB <- vector(mode="numeric", length=nmom)
       LA[3:nmom] <- NA; LB[1:2] <- NA
       for(i in 1:2) {
          LA[i] <- var(replicate(nsim,
                       lmoms(sample(x, n, replace=TRUE), nmom=nmom)$lambdas[i]))
       }
       for(i in 3:length(LB)) {
          LB[i] <- var(replicate(nsim,
                       lmoms(sample(x, n, replace=TRUE), nmom=nmom)$ratios[i]))
       }
       z <- list(lambdas=lmr$lambdas, ratios=lmr$ratios,
                 lambdavars=LA,       ratiovars=LB,
                 varcovar.lambdas=NA,
                 varcovar.lambdas.and.ratios=NA,
                 bootstrap.orderstatistics=NA,
                 varcovar.orderstatistics=NA,
                 inverse.varcovar.tau23=NA,
                 inverse.varcovar.tau34=NA,
                 inverse.varcovar.tau46=NA,
                 message="used simulation and not exact boot strap because of large sample size (n > bign)",
                 source="lmoms.bootbarvar")
        return(z)
    }

    L <- R <- V <- W <- vector(mode="numeric",length=nmom)
    weight.matrix <- matrix(nrow=nmom, ncol=n)
    if(nohatSIGMA) {
       hatsigma <- matrix(nrow=length(x), ncol=length(x))
    } else {
       hatsigma <- hatSIGMA(x) # a CPU hog
    }

    bootstrap.ostats <- hatMUk(x) #  the bootstrap estimates of the order statistics
    for(r in 1:nmom) {
       lmom.weights <- sapply(1:n, function(i) { Lcomoment.Wk(r, i, n)/n })
       lmom.weights <- as.vector(lmom.weights)
       weight.matrix[r,] <- lmom.weights
       bootbar <- sum(lmom.weights * bootstrap.ostats)
       bootvar <- crossprod(as.vector(crossprod(lmom.weights, hatsigma)), lmom.weights)
       L[r] <- bootbar
       V[r] <- bootvar
    }
    R[1] <- NA
    R[2] <- L[2]/L[1]
    R[3:nmom] <- L[3:nmom]/L[2]

    varcovarlam <- weight.matrix %*% hatsigma %*% t(weight.matrix)
    varcovartau <- matrix(nrow=nmom, ncol=nmom)
    #print(varcovarlam)
    #print(diag(varcovarlam))
    #print(V)

    for(r in 1:nmom) {
       for(s in 1:nmom) {
          #if(r < s) next
          if(r <= 2 && s <=2) {
             varcovartau[r,s] <- varcovarlam[r,s]
          } else if(r >= 3 && s <=2) {
             varcovartau[r,s] <- (varcovarlam[r,s] - lmr$ratios[r]*varcovarlam[2,s])/lmr$lambdas[2]
          } else if(r >= 3 && s >=3) {
             varcovartau[r,s] <- (varcovarlam[r,s] - lmr$ratios[r]*varcovarlam[2,s] -
                                                     lmr$ratios[s]*varcovarlam[2,r] +
                                  lmr$ratios[r]*lmr$ratios[s]*varcovarlam[2,2])/lmr$lambdas[2]^2
          }
       }
    }
    varcovartau[1,] <- varcovartau[,1]
    varcovartau[2,] <- varcovartau[,2]
    W <- diag(varcovartau)
    W[1:2] <- rep(NA, 2)
    if(covarinverse) {
       inv23 <- "L-CV and L-skew inversion not yet derived for same reason varcovar.lambdas.and.ratios[2]=NA"

       if(nmom < 4) {
          inv34 <- "too few L-moments requested in nmom (r < 4) for L-skew and L-kurtosis inversion"
       } else {
          tmp <- varcovartau[c(3,4), c(3,4)]
          inv34 <- chol2inv(chol(tmp))
       }
       if(nmom < 6) {
          inv46 <- "too few L-moments requested in nmom (r < 6) for L-kurtosis and Tau6 inversion"
       } else {
          tmp <- varcovartau[c(4,6), c(4,6)]
          inv46 <- chol2inv(chol(tmp))
       }

       z <- list(lambdas=L,    ratios=R,
                 lambdavars=V, ratiovars=W,
                 varcovar.lambdas=varcovarlam,
                 varcovar.lambdas.and.ratios=varcovartau,
                 bootstrap.orderstatistics=bootstrap.ostats,
                 varcovar.orderstatistics=hatsigma,
                 inverse.varcovar.tau23=inv23,
                 inverse.varcovar.tau34=inv34,
                 inverse.varcovar.tau46=inv46,
                 message="exact boot strap formula used",
                 source="lmoms.bootbarvar")
    } else {
       z <- list(lambdas=L,    ratios=R,
                 lambdavars=V, ratiovars=W,
                 varcovar.lambdas=varcovarlam,
                 varcovar.lambdas.and.ratios=varcovartau,
                 bootstrap.orderstatistics=bootstrap.ostats,
                 varcovar.orderstatistics=hatsigma,
                 inverse.varcovar.tau23=NA,
                 inverse.varcovar.tau34=NA,
                 inverse.varcovar.tau46=NA,
                 message="exact boot strap formula used",
                 source="lmoms.bootbarvar")
    }
    return(z)
}
