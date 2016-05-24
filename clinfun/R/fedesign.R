fe.power <- function(d, n1, n2, p1, alpha = 0.05) {
  nmx <- n1 + n2 + 1
  lgnmx <- lgamma(1:nmx)
  power <- 0
  zz1 <- .Fortran("ferej",
            as.integer(nmx),
            as.integer(n1),
            as.integer(n2),
            as.double(alpha),
            fcl = as.integer(matrix(0,nmx,2)),
            as.double(lgnmx),
            PACKAGE="clinfun")
  zz3 <- sapply(d,function(d,nmx,n1,n2,p1,fcl,lgnmx) {
    p2 <- p1+d
    zz2 <- .Fortran("fepow",
              as.integer(nmx),
              as.integer(n1),
              as.integer(n2),
              as.double(p1),
              as.double(p2),
              as.integer(fcl),
              as.double(lgnmx),
              power = as.double(double(1)),
              PACKAGE="clinfun")
    c(p2,zz2$power)
  },nmx,n1,n2,p1,zz1$fcl,lgnmx)
  zz3 <- t(zz3)
  zz3 <- as.data.frame(zz3)
  names(zz3) <- c("p2","power")
  zz3
}

fe.mdor <- function(ncase,ncontrol,pcontrol,alpha=0.05,power=0.8) {
  za <- qnorm(1-alpha/2)
  zb <- qnorm(power)
  cpc <- ncontrol/ncase
  A1 <- (za + zb)^2
  B1 <- 1 + 2 * cpc * pcontrol
  C1 <- 2 * cpc * pcontrol * (ncase * (1 - pcontrol) - A1 * pcontrol)
  if (C1 <= 0) { 
    warning("Small Sample Sizes and Large pcontrol") 
    smdor <- 1e8
  } 
  else {
    temp <- sqrt((A1 * B1)^2 + 2 * (1 + cpc) * A1 * C1)
    r1 <- 1 + (A1 * B1 - temp)/C1
    r2 <- 1 + (A1 * B1 + temp)/C1
    smdor <- max(r1, r2)
  }
  omdor <- matrix(0,3,2)
  omdor[1,1] <- smdor
  omdor[2,1] <- orCPS(ncase,ncontrol,pcontrol,alpha,power)
  n1 <- ncontrol
  n2 <- ncase
  nmx <- n1 + n2 + 1
  lgnmx <- lgamma(1:nmx)
  zz2 <- .Fortran("femdor",
            as.integer(nmx),
            as.integer(ncontrol),
            as.integer(ncase),
            as.double(pcontrol),
            as.double(alpha),
            as.double(power),
            as.integer(matrix(0,nmx,2)),
            as.double(lgnmx),
            omdor = as.double(omdor),
            PACKAGE="clinfun")
  xx <- matrix(zz2$omdor,3,2)
  dimnames(xx) <- list(c("Schlesselman","CPS","Fisher Exact"),c("Odds Ratio", "Exact Power"))
  xx
}

fe.ssize <- function(p1,p2,alpha=0.05,power=0.8,r=1,npm=5,mmax=1000) {
  m <- mCPS(p1,p2,alpha,power,r)
  if (m <= mmax) {
    ossiz <- matrix(0,2,3)
    ossiz[1, 1] <- m
    ossiz[1, 2] <- m*r
    nmx <- ceiling((m+npm)*(1+r)) + 2
    lgnmx <- lgamma(1:nmx)
    zz2 <- .Fortran("fessiz",
                    as.integer(nmx),
                    as.double(p1),
                    as.double(p2),
                    as.double(r),
                    as.double(alpha),
                    as.double(power),
                    as.integer(npm),
                    as.integer(matrix(0,nmx,2)),
                    as.double(lgnmx),
                    ossiz = as.double(ossiz),
                    PACKAGE="clinfun")
    xx <- matrix(zz2$ossiz,2,3)
    xx[1,1] <- ceiling(xx[1,1])
    xx[1,2] <- ceiling(xx[1,2])
    xx[2,1] <- ceiling(xx[2,1])
    xx[2,2] <- ceiling(xx[2,2])
    dimnames(xx) <- list(c("CPS","Fisher Exact"),c("Group 1", "Group 2", "Exact Power"))
  } else {
    message("\n Exact power is not computed for sample size greater than ", mmax, "\n")
    xx <- matrix(0,1,3)
    xx[1,1] <- ceiling(m)
    xx[1,2] <- ceiling(m*r)
    xx[1,3] <- power
    dimnames(xx) <- list(c("CPS"),c("Group 1", "Group 2", "Approx Power"))
  }
  xx  
}

mCPS <- function(p1,p2,alpha,power,r) {
  za <- qnorm(1-alpha/2)
  zb <- qnorm(power)
  p <- (p1+r*p2)/(1+r)
  m1 <- (za*sqrt((1+r)*p*(1-p)) + zb*sqrt(r*p1*(1-p1)+p2*(1-p2)))^2/(r*(p1-p2)^2)
  (m1*(1+sqrt(1+2*(r+1)/(m1*r*abs(p1-p2))))^2)/4
}

CPS.ssize <- function(p1,p2,alpha=0.05,power=0.8,r=1) {
  m <- mCPS(p1,p2,alpha=alpha,power=power,r=r)
  xx <- matrix(0, 1, 3)
  xx[1,1] <- m
  xx[1,2] <- m*r
  xx[1,1] <- ceiling(xx[1,1])
  xx[1,2] <- ceiling(xx[1,2])
  xx[1,3] <- power
  dimnames(xx) <- list(c("CPS"),c("Group 1", "Group 2", "Approx Power"))
  xx
}

orCPS <- function(ncase,ncontrol,pcontrol,alpha=0.05,power=0.8) {
  r <- ncase/ncontrol
  ff <- function(d, r, ncontrol, pcontrol, alpha, power) {
    za <- qnorm(1-alpha/2)
    zb <- qnorm(power)
    p1 <- pcontrol
    p2 <- pcontrol + d
    p <- (p1+r*p2)/(1+r)
    A <- (za*sqrt((1+r)*p*(1-p)) + zb*sqrt(r*p1*(1-p1)+p2*(1-p2)))^2
    m1 <- A/(r*(p1-p2)^2)
    (m1*(1+sqrt(1+2*(r+1)/(m1*r*abs(p1-p2))))^2)/4 - ncontrol
  }
  d0 <- uniroot(ff, lower=0.0001, upper=1-pcontrol-0.0001, r=r, ncontrol=ncontrol, pcontrol=pcontrol, alpha=alpha, power=power)$root
  (pcontrol+d0)*(1-pcontrol)/(pcontrol*(1-pcontrol-d0))
}

or2pcase <- function(pcontrol, OR) {OR*pcontrol/(1 + (OR -1)*pcontrol)}

# mdrr function added 2/18/2015
# estimate the minimum detectable response rate difference
# total sample size and marginal proportions are specified
# presp is one of the margins for example overall response rate of drug
# cprob is the other margin, say probability of a genetic aberration
# goal is to test if aberration is associated with response
mdrr <- function(n, cprob, presp, alpha=0.05, power=0.8, niter=15) {
    ncontrol <- round(n*cprob)
    ncase <- n-ncontrol
    out <- matrix(0,2,3)
    # overall rate presp = cprob*pcontrol + (1-cprob)*pcase
    # for numerical purposes lowest pcontrol is 1e-6 and highest pcase is 1-1e-6
    # pcontrol has to be as large as (presp - (1-cprob))/cprob
    pcontrol0 <- max(0.000001, (presp - 0.999999*(1-cprob))/cprob)
    pcase0 <- (presp - cprob*pcontrol0)/(1-cprob)
    pcontrol1 <- pcase1 <- presp
    opower0 <- fe.power(pcase0-pcontrol0, ncontrol, ncase, pcontrol0, alpha)$power
    opower1 <- fe.power(pcase1-pcontrol1, ncontrol, ncase, pcontrol1, alpha)$power
    if (opower0 > power) {
        # power at mid-points
        for(i in 1:niter) {
            pcontrol <- (pcontrol0+pcontrol1)/2
            pcase <- (presp - cprob*pcontrol)/(1-cprob)
            opower <- fe.power(pcase-pcontrol, ncontrol, ncase, pcontrol, alpha)$power
            if (opower < power) {
                pcontrol1 <- pcontrol
                pcase1 <- pcase
                opower1 <- opower
            } else {
                pcontrol0 <- pcontrol
                pcase0 <- pcase
                opower0 <- opower
            }
        }
    } else {
        pcontrol <- pcontrol0
        pcase <- pcase0
        opower <- opower0        
        warning("required power cannot be achieved")
    }
    out[1,] <- c(pcontrol, pcase, opower) # pclass > presp
    # repeat now to get pclass < presp
    pcontrol0 <- min(0.999999, (presp - 0.000001*(1-cprob))/cprob)
    pcase0 <- (presp - cprob*pcontrol0)/(1-cprob)
    pcontrol1 <- pcase1 <- presp
    opower0 <- fe.power(pcase0-pcontrol0, ncontrol, ncase, pcontrol0, alpha)$power
    opower1 <- fe.power(pcase1-pcontrol1, ncontrol, ncase, pcontrol1, alpha)$power
    if (opower0 > power) {
        # power at mid-points
        for(i in 1:niter) {
            pcontrol <- (pcontrol0+pcontrol1)/2
            pcase <- (presp - cprob*pcontrol)/(1-cprob)
            opower <- fe.power(pcase-pcontrol, ncontrol, ncase, pcontrol, alpha)$power
            if (opower < power) {
                pcontrol1 <- pcontrol
                pcase1 <- pcase
                opower1 <- opower
            } else {
                pcontrol0 <- pcontrol
                pcase0 <- pcase
                opower0 <- opower
            }
        }
    } else {
        pcontrol <- pcontrol0
        pcase <- pcase0
        opower <- opower0        
        warning("required power cannot be achieved")
    }
    out[2,] <- c(pcontrol, pcase, opower)
    colnames(out) <- c("P(resp|+)","P(resp|-)", "power")
    out
}
