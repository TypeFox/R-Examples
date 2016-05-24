HWData <- function (nm = 100, n = rep(100, nm), f = rep(0, nm), p = NULL, 
                    conditional = FALSE, exactequilibrium = FALSE, 
                    pdist = "runif", x.linked = FALSE,
                    nA = NULL, n.males=rep(round(0.5*n),nm),...) 
{
    if (length(p) == 1) 
        p <- rep(p, nm)
    if (length(n) == 1) 
        n <- rep(n, nm)
    if (length(f) == 1) 
        f <- rep(f, nm)
     if (length(nA) == 1) 
        nA <- rep(nA, nm)
    if (length(n.males) == 1) 
        n.males <- rep(n.males, nm)
    nAll <- 2*n
    # setup the allele frequencies or allele counts
    if(any(is.null(nA))) { # look at p if nA not given
      if(any(is.null(p)))  # take at random of p not given
        if (pdist == "runif") 
          p <- runif(nm, ...)
        else {
          if (pdist == "rbeta") 
            p <- rbeta(nm, ...)
          else stop("unknown value for pdist")
        }        
      if(!x.linked) nA <- round(p * nAll, digits = 0) else {
        nt <- n.males + 2*(n-n.males)
        nA <- round(p*nt, digits = 0) 
        }# set the allele counts if unspecified
    } else { # if not null nA is specified and we compute allele frequencies.
      if(!x.linked) p <- nA/(2*n) else {
        nf <- n - n.males
        nt <- 2*nf + n.males
        p <- nA/nt
      }
    }
    if(!x.linked) { # autosomal marker
      X <- matrix(0,nrow=nm,ncol=3)
      colnames(X) <- c("AA","AB","BB")
      if(!conditional) { # multinomial sampling
        if(!exactequilibrium) {
          for(i in 1:nm) {
            X[i,] <- t(rmultinom(1, size = n[i], 
                            prob = c(p[i]^2 + f[i] * p[i] * (1 - p[i]), 
                                    (1 - f[i]) * 2 * p[i] * (1 - p[i]), 
                                    (1 - p[i])^2 + f[i] * p[i] * (1 - p[i]))))
          }
        } else { # simulate data in exact equilibrium
          for(i in 1:nm) {
            X[i,] <- c(n[i] * p[i]^2 + n[i] * f[i] * p[i] * (1 - p[i]), 
                   n[i] * (1 - f[i]) * 2 * p[i] * (1 - p[i]), 
                   n[i] * (1 - p[i])^2 + n[i] * f[i] * p[i] * (1 - p[i]))
          }
        }
      } else { # conditional Levene-Haldane sampling 
        nB <- nAll - nA
        for (i in 1:nm) {
          Pop <- c(rep(1, nA[i]), rep(0, nB[i]))
          sam <- matrix(sample(Pop, nAll), ncol = 2, byrow = TRUE)
          status <- apply(sam, 1, sum)
          nAA <- sum(status == 2)
          nAB <- sum(status == 1)
          nBB <- sum(status == 0)
          if ((nAA + nAB + nBB) != n[i]) 
            stop("HWData: error")
          X[i,] <- c(nAA, nAB, nBB)
        }
      }
    } 
    else {# if (!x.linked) (marker on the X chromosome)
      theta <- n.males/n
      nf <- n - n.males
      nt <- 2*nf + n.males
      nB <- nt - nA
      if(!all(nA <= nB)) warning("Major allele specified instead of Minor allele")
      pA <- nA/nt
      prob <- cbind(A=theta*pA,theta*(1-pA),(1-theta)*pA^2,(1-theta)*2*pA*(1-pA),(1-theta)*(1-pA)^2)
      X <- matrix(0,nrow=nm,ncol=5)
      colnames(X) <- c("A","B","AA","AB","BB")
      if(!conditional) {
        for(i in 1:nm) {
          row <- t(rmultinom(1,size=n[i],prob[i,])) 
          X[i,] <- row
        }
      } else {
          for(i in 1:nm) {
          Pop <- c(rep(1,nA[i]),rep(0,nt[i]-nA[i]))
          Sam <- sample(Pop,nt[i])
          Males <- Sam[1:n.males[i]]
          mA <- sum(Males==1)
          mB <- sum(Males==0)
          Females <- matrix(Sam[(n.males[i]+1):nt[i]],ncol=2,byrow=TRUE)      
          status <- apply(Females, 1, sum)
          fAA <- sum(status == 2)
          fAB <- sum(status == 1)
          fBB <- sum(status == 0)
          if ((fAA + fAB + fBB) != nf[i]) stop("HWData: error")
          X[i,] <- c(mA, mB, fAA, fAB, fBB)     
      }
    }
  } # if !x.linked
  if(nm==1) { # return a named vector if only one marker is asked for.
      c.names <- colnames(X)
      X <- as.vector(X[1,])
      names(X) <- c.names
  }
  return(X)
}




