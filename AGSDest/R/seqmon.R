seqmon <- function (a, b, t, int) {
    ones <- function(a, b) {
        matrix(1, nrow=a, ncol=b)
        ##array(rep(1, a * b), c(a, b))
    }
        
    d <- (b - a)/int
    m <- length(a)
    pU = ones(m, 1)
    pL = ones(m, 1)
    sq2pi <- sqrt(2 * pi)
    H <- 1:int[1]
    E <- ones(1, int[1])
    
    xo <- a[1] + ((1:int[1]) - 0.5 * E) * d[1]
    
    pU[1] <- pnorm(-(sqrt(t[1]) * b[1])/sqrt(t[1]))
    
    M <- t((d[1]/sq2pi) * exp(-(sqrt(t[1]) * xo)^2/(2 * t[1])))
    
    pL[1] <- pnorm(sqrt(t[1]) * a[1]/sqrt(t[1]))
    
    for(k in 2:m) {
        VU <- pnorm(-(sqrt(t[k]) * b[k] * E - sqrt(t[k - 1]) * xo)/sqrt(t[k] - t[k - 1]))
        VL <- pnorm((sqrt(t[k]) * a[k] * E - sqrt(t[k - 1]) * xo)/sqrt(t[k] - t[k - 1]))
        pL[k] <- pL[k - 1] + VL %*% M
        pU[k] <- pU[k - 1] + VU %*% M
        x <- a[k] + ((1:int[k]) - 0.5 * ones(1, int[k])) * d[k]
        
        if(k != m) {
            M <- (d[k] * sqrt(t[k])/(sq2pi * sqrt(t[k] - t[k - 1]))) * exp(-(sqrt(t[k]) * (t(x) %*% ones(1, int[k - 1])) - sqrt(t[k - 1]) * (ones(int[k], 1) %*% xo))^2/(2 * (t[k] - t[k - 1]))) %*% M
            xo <- x
        }
    }
    
    c(pL, pU)
}

A <- function(h, k=length(pT$t), pT, iD) {
    if(k-iD$T <= 1) 0
    else {
        if(k-iD$T == 2) 1-pnorm(pbounds(h, pT, iD)[1])
        else {
            seqmon(a=pT$a[(iD$T+1):(k-1)],
                   b=pbounds(h,pT,iD)[1:(k-1-iD$T)],
                   t=pT$t[(iD$T+1):(k-1)]*pT$Imax-pT$t[iD$T]*pT$Imax,
                   int=500*array(c(1),k-1-iD$T))[2*(k-1-iD$T)]
        }
    }
}

comp.alab <- function(GSD, pprec=1e-05) {
    K <- length(GSD$a)
    alab <- numeric(K+1)
    for(k1 in 1:(K-1)) {
        if(k1 == 1) alab[1] <-(GSD$b[1] - qnorm(1-GSD$al)) / sqrt(GSD$t[1]*GSD$Imax)
        else {
            alab[k1] <- uniroot(function(x) seqmon(a=GSD$a[1:k1],
                                                   b=GSD$b[1:k1] - x*sqrt(GSD$t[1:k1]*GSD$Imax),
                                                   t=GSD$t[1:k1]*GSD$Imax,
                                                   int=500*rep.int(1, k1))[2*k1]-GSD$al,
                                c(0, alab[k1-1]))$root
        }
    }
    alab[K] <- 0
    ##lower critical boundary for delta's where the exceeding probability 
    ##(estimated by Bonferroni) is smaller than pprec
    alab[K+1] <-min((GSD$b[1:(K-1)] - qnorm(1-pprec/(1:(K-1)))) / sqrt(GSD$t[1:(K-1)]*GSD$Imax))
    alab[K+2] <-pprec;
    alab
}


j.alab <- function(h, GSD=NULL) {
    if(is.null(GSD)) stop("missing input: specify GSD") ## FIXME
    if(is.null(GSD$alab)) GSD$alab <- comp.alab(GSD=GSD)
    length(GSD$alab[1:(length(GSD$alab)-3)]) + 1 - sum(ifelse(h >= GSD$alab[1:(length(GSD$alab)-3)], 1, 0))
}

bh <- function(h, pT, xl=NULL, xu=NULL) {
    if(h == 0) bh <- pT$b[length(pT$t)]
    else {
        if(is.null(pT$alab)) pT$alab <- comp.alab(GSD=pT)
        if(h %in% pT$alab[1:(length(pT$t)-1)]) Inf
        else {
            k  <- j.alab(h, GSD=pT)
            if(is.null(xl)) xl <- qnorm(1-pT$al) + h * sqrt(pT$t[k]*pT$Imax)
            if((k == 1) | (h <= pT$alab[length(pT$t) + 1])) xl
            else {
                if(k < length(pT$t)) xl <-max(xl, pT$b[k])
                if(is.null(xu)) {
                    al0 <- ifelse(k == 2, 1-pnorm(pT$b[(k-1)]- h * sqrt(pT$t[(k-1)] * pT$Imax)),
                                  seqmon(a=pT$a[1:(k-1)], b=pT$b[1:(k-1)]-h*sqrt(pT$t[1:(k-1)]*pT$Imax), t=pT$t[1:(k-1)] * pT$Imax,
                                         int=500*rep.int(1, k-1))[2*(k-1)]);
                    al1 <- (pT$al - al0) / (1-al0)
                    xu  <- qnorm(1-al1) + h * sqrt(pT$t[k] * pT$Imax)
                }
                uniroot(f=function(x) seqmon(a=pT$a[1:k], b=c(pT$b[1:(k-1)],x)-h*sqrt(pT$t[1:k]*pT$Imax), t=pT$t[1:k]*pT$Imax,
                            int=500*rep.int(1, k))[2*k]-pT$al, interval=c(xl,xu))$root
              }
        }
    }
}

comp.als <- function(GSD) { 
    K <- length(GSD$t)
    if(K == 1) 1-pnorm(GSD$b)
    else seqmon(a=GSD$a, b=GSD$b, t=GSD$t*GSD$Imax, int=500*rep.int(1, K))[(K+1):(2*K)]
}

compBounds <- function(t=1:length(t2)/length(t2), t2, iuse, asf = NULL, alpha, phi, ztrun) {
#   if(whatSpendingFunctionIsUsed==1 || whatSpendingFunctionIsUsed==3){
   K <- length(t)
   n <- length(t)
   if(iuse == 3) {
      bmin <- qnorm(1-alpha)
      bmax <- qnorm(1-(alpha/n))/min(1,n^(phi-0.5))
      if(n > 1) {
          b <- uniroot(f=function(x)seqmon(a=rep(-ztrun,n),
                           b=x*(1:n)^(phi-0.5),
                           t=t2,
                           int=500*array(1,n))[2*n]-alpha
                      ,interval=c(bmin,bmax))$root;
          b*(1:n)^(phi-0.5)
      } else qnorm(1-alpha)
  } else {
      if(iuse == 4) {
          if(phi == 0) stop("phi must be unequal 0")
          else {
              if(n > 1) {
                  levsp <- alpha*(1-exp(-phi*t2/t2[n]))/(1-exp(-phi))
                  b <- numeric(n)
                  b[1] <- qnorm(1-levsp[1])
                  for(i in 2:n) {
                      bmin <- qnorm(1-levsp[i])
                      bmax <- qnorm(1-levsp[i]+levsp[i-1])
                      if(abs(bmin-bmax)<=10^(-4)) b[i] <- bmin
                      else {
                          b[i] <- uniroot(function(x) seqmon(a=rep(-ztrun,i),
                                                             b=c(b[1:(i-1)],x),
                                                             t=t2[1:i],
                                                             int=500*array(1,i))[2*i]-levsp[i],
                                          interval=c(bmin,bmax))$root
                      }
                  }
                  b
              } else qnorm(1-alpha)
          }
      } else ldbounds::bounds(t=1:K/K, t2 = t2, iuse = iuse, asf = NULL,alpha = alpha, phi = phi,ztrun = 8)$upper.bounds
  }
}
