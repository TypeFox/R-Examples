dskellam <- function(x, lambda1, lambda2=lambda1, log=FALSE){
 # density (PMF) of Skellam distriubition (difference of Poissons)
    if (missing(x)|missing(lambda1)) stop("first 2 arguments are required")
    lambdas <- c(lambda1,lambda2)
    oops <- !(is.finite(lambdas)&(lambdas>=0))
    if(any(oops)) {
        warning("NaNs produced")
        lambdas[oops] <- NaN
        lambda1 <- lambdas[1:length(lambda1)]
        lambda2 <- lambdas[(length(lambda1)+1):length(lambdas)]
    }
   # make all args the same length (for subsetting)
    lens <- c(length(x),length(lambda1),length(lambda2))
    len <- max(lens)
    if(len>min(lens)) {
        if (all(len/lens!=len%/%lens)) warning("longer object length is not a multiple of shorter object length", domain=NA)
        x <- rep(x,length.out=len)
        lambda1 <- rep(lambda1,length.out=len)
        lambda2 <- rep(lambda2,length.out=len)
    }
   # warn of non-integer x values (since support of PMF is integers)
    nonint <- x!=trunc(x)
    if (any(nonint)) {
        xreal <- x[nonint]
        for (i in 1:length(xreal)) warning(paste("non-integer x =",xreal[i]))
    }
    x <- trunc(x)
    ret <- rep(NaN,length.out=len)
   # handle a zero lambda separately (avoids division by zero & accuracy issues for large values of lambda or x)
    ret[lambda1==0] <- dpois(-x[lambda1==0],lambda2[lambda1==0],log=log)
    ret[lambda2==0] <- dpois( x[lambda2==0],lambda1[lambda2==0],log=log)    # corrects bug in VGAM 0.7-9
   # non-zero lambdas
    nz <- is.finite(lambda1)&is.finite(lambda2)&(lambda1>0)&(lambda2>0)
    L1 <- lambda1[nz];
    L2 <- lambda2[nz];
    sqL12 <- sqrt(L1*L2)
    xx <- x[nz]
    if (log[1]) {       # use abs(x) in besselI for improved accuracy (prior to 2.9) and S-PLUS compatibility
       # log(besselI(y,nu)) == y+log(besselI(y,nu,TRUE))
        ret[nz] <- log(besselI(2*sqL12, abs(xx),TRUE))+2*sqL12-L1-L2+xx/2*log(L1/L2)
    } else {
       # besselI(y,nu); exp(y)*besselI(y,nu,TRUE)
#       ret[nz] <- besselI(2*sqL12, abs(xx),TRUE)*(L1/L2)^(xx/2)*exp(2*sqL12-L1-L2)
        ret[nz] <- besselI(2*sqL12, abs(xx),TRUE)*exp(2*sqL12-L1-L2+xx/2*log(L1/L2))
    }
    chk <- nz & (!is.finite(ret) | (!log)&(ret<1e-308)) # use saddlepoint approximation to detect possible over/underflow
    if (length(chk[chk])>0) {
        L1 <- lambda1[chk];
        L2 <- lambda2[chk];
            sqL12 <- sqrt(L1*L2)
        xx <- x[chk]
        s <- log(0.5*(xx+sqrt(xx^2+4*L1*L2))/L1)# the saddlepoint
        K <- L1*(exp(s)-1)+L2*(exp(-s)-1)       # CGF(s)
        K2 <- L1*exp(s)+L2*exp(-s)              # CGF''(s)
        spd <- exp(K-x*s)/sqrt(2*pi*K2)         # saddlepoint density
        usp <- (spd>1e-308)&is.finite(spd)      # don't trust the existing result
        if (length(usp[usp])>0) {   # add another term to the saddlepoint approximation
            su <- s[usp]
            K2u <- K2[usp]
            c <- (1-((L1[usp]*exp(su)-L2[usp]*exp(-su))/K2u)^2*5/3)/K2u*0.125+1
            ret[chk][usp] <- exp(K[usp]-x[usp]*su)/sqrt(2*pi*K2u)*(1+c)*0.5
        }
    }
    ret
}
