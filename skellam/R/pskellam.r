pskellam <- function(q, lambda1, lambda2=lambda1, lower.tail=TRUE, log.p=FALSE){
 # CDF of Skellam distriubition (difference of Poissons)
    if (missing(q)|missing(lambda1)) stop("first 2 arguments are required")
    lambdas <- c(lambda1,lambda2)
    oops <- !(is.finite(lambdas)&(lambdas>=0))
    if(any(oops)) {
        warning("NaNs produced")
        lambdas[oops] <- NaN
        lambda1 <- lambdas[1:length(lambda1)]
        lambda2 <- lambdas[(length(lambda1)+1):length(lambdas)]
    }
   # CDF is a step function, so convert to integer values without warning
    x <- floor(q)
   # make all args the same length (for subsetting)
    lens <- c(length(x),length(lambda1),length(lambda2))
    len <- max(lens)
    if(len>min(lens)) {
        if (all(len/lens!=len%/%lens)) warning("longer object length is not a multiple of shorter object length", domain=NA)
        x <- rep(x,length.out=len)
        lambda1 <- rep(lambda1,length.out=len)
        lambda2 <- rep(lambda2,length.out=len)
    }
   # different formulas for negative & nonnegative x (zero lambda is OK)
    neg <- (x< 0)&(!is.nan(lambda1))&(!is.nan(lambda2))
    pos <- (x>=0)&(!is.nan(lambda1))&(!is.nan(lambda2))
    ret <- rep(NaN,length.out=len)
    if (lower.tail[1]) {
        ret[neg] <- pchisq(2*lambda2[neg],-2*x[neg],2*lambda1[neg],log.p=log.p)
        ret[pos] <- pchisq(2*lambda1[pos],2*(x[pos]+1),2*lambda2[pos],lower.tail=FALSE,log.p=log.p)     # S-PLUS does not have a lower.tail argument
    } else {
        ret[neg] <- pchisq(2*lambda2[neg],-2*x[neg],2*lambda1[neg],lower.tail=FALSE,log.p=log.p)        # S-PLUS does not have a lower.tail argument
        ret[pos] <- pchisq(2*lambda1[pos],2*(x[pos]+1),2*lambda2[pos],log.p=log.p)
    }
    chk <- (neg|pos) & (!is.finite(ret) | (!log.p)&(ret<1e-308))    # use saddlepoint approximation if outside the working range of pchisq
    if (length(chk[chk])>0) ret[chk] <- pskellam.sp(x[chk],lambda1[chk],lambda2[chk],lower.tail,log.p)
    ret
}
