rskellam <- function(n, lambda1, lambda2=lambda1){
 # Skellam random variables
    if (missing(n)|missing(lambda1)) stop("first 2 arguments are required")
    if (length(n)>1) n <- length(n)
    lambda1 <- rep(lambda1,length.out=n)
    lambda2 <- rep(lambda2,length.out=n)
    oops <- !(is.finite(lambda1)&(lambda1>=0)&is.finite(lambda2)&(lambda2>=0))
    if(any(oops)) warning("NaNs produced")
    ret <- rep(NaN,length.out=n)
    n <- n-sum(oops)
    ret[!oops] <- rpois(n,lambda1[!oops])-rpois(n,lambda2[!oops])
    ret
}
