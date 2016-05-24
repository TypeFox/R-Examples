#binomial coefficient
binom.coef=function(n,m) { prod((n-m+1):n)/prod(1:m) }

expit=function (x) {1/(1+exp(-x))}
logit=function (x) {log(x/(1-x))}


as.binary <- function(n, base=2 , r=FALSE) {
   out <- NULL
   while(n > 0) {
     if(r) {
       out <- c(out , n%%base)
     } else {
       out <- c(n%%base , out)
     }   
     n <- n %/% base
   }
   return(out)
}



##-- Stirling numbers of the 2nd kind
##-- (Abramowitz/Stegun: 24,1,4 (p. 824-5 ; Table 24.4, p.835)

##> S^{(m)}_n = number of ways of partitioning a set of $n$ elements into $m$
##> non-empty subsets

Stirling2 <- function(n,m)
{
    ## Purpose:  Stirling Numbers of the 2-nd kind
    ##      S^{(m)}_n = number of ways of partitioning a set of
    ##                      $n$ elements into $m$ non-empty subsets
    ## Author: Martin Maechler, Date:  May 28 1992, 23:42
    ## ----------------------------------------------------------------
    ## Abramowitz/Stegun: 24,1,4 (p. 824-5 ; Table 24.4, p.835)
    ## Closed Form : p.824 "C."
    ## ----------------------------------------------------------------
    ## maechler@_stat.math.ethz.ch
    
    if (0 > m || m > n) stop("'m' must be in 0..n !")     
    k <- 0:m
    sig <- rep(c(1,-1)*(-1)^m, length= m+1)
    # 1 for m=0; -1 1 (m=1)     
    ## The following gives rounding errors for (25,5) :     
    ## r <- sum( sig * k^n /(gamma(k+1)*gamma(m+1-k)) )     
    ga <- gamma(k+1)
    round(sum( sig * k^n /(ga * rev(ga)))) 
}

logSumExp=function (logx){
    logMeanExp(logx, 1)
}
logSumExpFor2=function (logx, logy){
    c=max(logx, logy)
    dif=abs(logx-logy)
    if (dif>300) return (c)
    else {
        log(sum(exp(logx-c), exp(logy-c)))+c
    }
}
# log( exp(logx1)-exp(logx2) )
logDiffExp=function (logx1,logx2){
    c=logx1
    c+ log(1-exp(logx2-logx1))
}

logMeanExp=function (logx,B=NULL){
# mean function for small numbers
# logx is a vector of large negative values
# return log (sum(exp(logx))/B)
# calculate log of the mean of a vector which contains 0 and very small real numbers (logged)
# return log of the mean
    if (is.null(B)) B=length(logx)
    c=max(logx)
    log(sum(exp(logx-c))/B)+c
}
# logMeanExp(log(1:5), 5) # test, should return log(3)

logDiffExp=function (logx1, logx2){
# diff function for small numbers
# logx1 and logx2 are typically large negative values, logx1>logx2
# return log (exp(logx1)-exp(logx2))
    if (logx1<logx2) {cat("\nlWarning [logDiffExp]: first argument smaller than second, return NaN.\n\n"); return (NaN);}
    c=max(logx1, logx2)
    log (exp(logx1-c)-exp(logx2-c))+c
}
# logDiffExp(log(2), log(1)) # test, should return 0

# from combinat package
permn=function (x, fun = NULL, ...) 
{
    if (is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == 
        x) 
        x <- seq(x)
    n <- length(x)
    nofun <- is.null(fun)
    out <- vector("list", gamma(n + 1))
    p <- ip <- seqn <- 1:n
    d <- rep(-1, n)
    d[1] <- 0
    m <- n + 1
    p <- c(m, p, m)
    i <- 1
    use <- -c(1, n + 2)
    while (m != 1) {
        out[[i]] <- if (nofun) 
            x[p[use]]
        else fun(x[p[use]], ...)
        i <- i + 1
        m <- n
        chk <- (p[ip + d + 1] > seqn)
        m <- max(seqn[!chk])
        if (m < n) 
            d[(m + 1):n] <- -d[(m + 1):n]
        index1 <- ip[m] + 1
        index2 <- p[index1] <- p[index1 + d[m]]
        p[index1 + d[m]] <- m
        tmp <- ip[index2]
        ip[index2] <- ip[m]
        ip[m] <- tmp
    }
    out
}
