#
#  EM algorithm for Mixture of Unrestricted Multivariate Skew t-distributioins
#  Package: EMMIX-uskew
#  Version: 0.11-5
#
#  Code by S. X. Lee
#  Updated on 11 Mar, 2013
#
# Genz, A. and Bretz, F. (2002) Comparison of Methods for the Numerical
#   Computation of Multivariate t Probabilities
#

################################################################################
#  SECTION 5
#                 Multivariate t Distribution Function
#                  for non-intger degrees of freedom
#
################################################################################
#
# Implementation of a method described in:
#     "Comparison of Methods for the Numerical Computation of Multivariate
#     t Probabilities", J. of Comput. Graph. Stat., 11(2002), pp. 950-971,
#        by Alan Genz and Frank Bretz.
#

mtcdf <- function(dat, mu = rep(0, length(dat)), sigma = diag(length(dat)), dof = Inf, ...){
      m <- 1e4*length(dat)
      p <- length(dat)
      cn <- diag(p)
      a <- rep(-Inf, p)
      return(qscmvtv(m, dof, sigma, a, cn, dat-mu))
}

mtcdf2 <- function (dat, mu = rep(0, length(dat)), sigma = diag(length(dat)), dof = Inf, ...) {
    T1 <- pmt(dat, mu, sigma, floor(dof), method=1, ...)
    T2 <- pmt(dat, mu, sigma, ceiling(dof), method=1, ...)
    rem <- dof %% 1
    return((rem*T2 + (1-rem)*T1))
}

qscmvtv <- function(m, nu, sigma, a, cn, b) {
    tmp <- chlrst(sigma, a, cn, b);
    as <- tmp[[1]]; ch <- tmp[[2]]; bs <- tmp[[3]]; clg <- tmp[[4]]; n <- tmp[[5]]
    p <- e <-  0;
    ns <- 8;
    nv <- trunc(max(m/(2*ns),1))                     
    if (nu > 0) nm <- n   else  nm <- n-1            
    on <- sp <- sm <- rep(1, nv)                     
    q <- 2^((1:nm)/(nm+1))                           
    for(i in 1:ns) {
        xx <- abs(2*((q%*%t(1:nv) + runif(nm)%*%t(on)) %% 1) -1)    
        if (nu > 0) {
            sp <- sqrt(2*qgamma(xx[n,], nu/2) / nu)
            sm <- sqrt(2*qgamma(1-xx[n,], nu/2) / nu)
        }
        vp <- mvtdnv(n, as, ch, bs, clg, xx, sp, nv)
        vp <- (mvtdnv( n, as, ch, bs, clg, 1-xx, sm, nv ) + vp )/2
        d <- (mean(vp) - p)/i                        
        p <- p + d
        if(abs(d) > 0)  e <- abs(d)*sqrt(1 + (e/d)^2 * (i-2)/i)
        else if(i>1) e <- e*sqrt((i-2)/i)
    }
    e <- 3*e
    return(c(p, e))
}

mvtdnv <- function(n, a, ch, b, clg, x, r, nv) {
    cc <- phi(a[1]*r)
    dc <- phi(b[1]*r) - cc
    p <- dc
    y <- matrix(0, n-1, nv)
    li <- 2; lf <- 1
    if(n>=2) {
    for(i in 2:n) {
        y[i-1,] <- phiinv(cc + x[i-1,]*dc)
        lf <- lf + clg[i]
        if(lf < li)   {cc <- 0;  dc <- 1}
        else {
            s <- matrix(ch[li:lf, 1:(i-1)], length(li:lf)) %*% y[1:(i-1),]
            ai <- matrix(a[li:lf])%*%r-s
            ai <- pmax(ai,-9)
            bi <- matrix(b[li:lf])%*%r - s;
            bi <- pmin(bi, 9)
            bi <- pmax(ai, bi)
            cc <- phi(ai)
            dc <- phi(bi) - cc
            p <- p*dc
        }
        li <- li + clg[i]
    }}
    return(p)
}


chlrst <- function(r, a, cn, b) {
    ep <- 1e-10   
    m <- nrow(cn);  n <- ncol(cn)
    ch <- cn
    np <- 0
    ap <- a; bp <- b
    y <- rep(0,n)
    sqtp <- sqrt(2*pi)
    cc <- r
    d <- sqrt((diag(cc)>0)*diag(cc) + (diag(cc)<0)*0)
    for(i in 1:n) {
        di <- d[i]
        if(di > 0)  {
            cc[,i] <- cc[,i]/di
            cc[i,] <- cc[i,]/di
            ch[,i] <- ch[,i]*di
        }
    }    
    clg <- rep(0, n)
    for(i in 1:n){
        epi <- ep*i
        j <- i
        if(i < n) {for(l in (i+1):n) {if(cc[l,l] > cc[j,j]) j <- l}}
        if(j > i) {
            tmp <- cc[i,i]; cc[i,i] <- cc[j,j]; cc[j,j] <- tmp;
            if(i > 1) {tmp <- cc[i,1:(i-1)]; cc[i,1:(i-1)] <- cc[j,1:(i-1)]; cc[j,1:(i-1)] <- tmp;}
            if((j-i)>2) {tmp <- cc[(i+1):(j-1),i]; cc[(i+1):(j-1),i] <- t(cc[j,(i+1):(j-1)]); cc[j,(i+1):(j-1)] <- t(tmp);}
            if(j < n) {tmp <- cc[(j+1):n,i]; cc[(j+1):n,i] <- cc[(j+1):n,j]; cc[(j+1):n,j] <- tmp;}
            tmp <- ch[,i]; ch[,i] <- ch[,j]; ch[,j] <- tmp;
        }
        if(cc[i,i]<epi) break
        cvd <- sqrt(cc[i,i]); cc[i,i] <- cvd
        if(i < n) {
            for(l in (i+1):n) {
                cc[l,i] <- cc[l,i]/cvd
                cc[l,(i+1):l] <- cc[l,(i+1):l] - cc[l,i]*t(cc[(i+1):l,i])
            }
        }
        if(i==n) ch[,i] <- ch[,i:n] * cc[i:n,i]  
        else ch[,i] <- ch[,i:n] %*% cc[i:n,i]
        np <- np + 1
    }
    if(min(np-1,m) >= 1){ 
    for(i in 1:min(np-1,m)) {
        epi <- epi*i
        vm <- 1
        lm <- i
        for(l in i:m) {
            v <- ch[l, 1:np]
            s <- {if(i>1) as.numeric(v[1:(i-1)]%*%y[1:(i-1)]) else 0}     
            ss <- max(sqrt(sum(v[i:np]^2)),epi)                           
            al <- (ap[l]-s)/ss; bl <- (bp[l]-s)/ss                        
            dna <- dsa <- dnb <- 0; dsb <- 1
            if(al > (-9)) {dna <- exp(-al*al/2)/sqtp; dsa <- phi(al)} else al <- bl
            if(bl < 9) {dnb <- exp(-bl*bl/2)/sqtp; dsb <- phi(bl)} else bl <- al
            if((dsb-dsa)>epi) {
                ds <- dsb - dsa
                mn <- (dna-dnb)/ds
                vr <- 1 + (al*dna -bl*dnb)/ds -mn^2
            }else{ mn <- (al+bl)/2; vr <- 0}
            if(vr <= vm) {lm <- l; vm <- vr; y[i] <- mn}
        }
        v <- ch[lm, 1:np]
        if(lm > i){
            ch[lm, 1:np] <- ch[i, 1:np]
            ch[i,1:np] <- v
            tl <- ap[i]; ap[i] <- ap[lm]; ap[lm] <- tl
            tl <- bp[i]; bp[i] <- bp[lm]; bp[lm] <- tl
        }
        if(i < np) ch[i, (i+1):np] <- 0
        if(i < np) ss <- sum(v[(i+1):np]^2) else ss <- 0
        if(ss > epi) {
            ss <- sqrt(ss+v[i]^2)
            if(v[i] < 0)  ss <- -ss
            ch[i,i] <- -ss
            v[i] <- v[i] + ss
            vt <- matrix((v[i:np])/(ss*v[i]),length(i:np))
            ch[(i+1):m, i:np] <- ch[(i+1):m, i:np] - ch[(i+1):m, i:np]%*%vt%*%v[i:np]
        }
    }}
    clm <- rep(0,m)
    for(i in 1:m) {
        v <- ch[i, 1:np]
        clm[i] <- min(i, np)
        jm <- 1;
        for(j in 1:clm[i])  {if(abs(v[j])> (ep*j)) jm <- j}
        if((jm+1)<=np) v[(jm+1):np] <- 0        
        clg[jm] <- clg[jm] + 1
        at <- ap[i]; bt <- bp[i]; j <- i
        if(i > 1) {
            for(l in (i-1):1) {
                if (jm >= clm[l]) break
                ch[(l+1), 1:np] <- ch[l, 1:np];  j <- 1
                ap[l+1] <- ap[l];  bp[l+1] <- bp[l];  clm[l+1] <- clm[l]
            }
        }       
        clm[j] <- jm; vjm <- v[jm]
        ch[j,1:np] <- v/vjm
        ap[j] <- at/vjm;  bp[j] <- bt/vjm
        if(vjm<0) {tl <- ap[j]; ap[j] <- bp[j]; bp[j] <- tl}
    }
    j <-0; for(i in 1:np) {if(clg[i]>0) j <- i}
    np <- j
    if(clg[1]>1) {
        ap[1] <- max(ap[1:clg[1]])
        bp[1] <- max(ap[1], min(bp[1:clg[1]]))
        ap[2:(m-clg[1]+1)] <- ap[(clg[1]+1):m]
        bp[2:(m-clg[1]+1)] <- bp[(clg[1]+1):m]
        clg[1] <- 1
    }
    return(list(ap,ch, bp, clg, np))
}

phi <- function(z) { 
    erfc <- function (x) 2 * pnorm(x * sqrt(2), lower.tail = FALSE)
    return(erfc(-z/sqrt(2))/2)
}

phiinv <- function(z) { 
    erfcinv <- function (x) qnorm(x/2, lower.tail = FALSE)/sqrt(2)
    return(-sqrt(2)*erfcinv(2*z))
}



