# $Id: ptte.data.R,v 1.2 2005/09/26 02:29:00 andrewr Exp $

"ptte.data" <-
function(x, alpha=0.05, Epsilon=0.25) {
    mean <- mean(x, na.rm=TRUE); std <- sd(x, na.rm=TRUE)
    n <- length(x) - sum(is.na(x))
    df1 <- 1; df2 <- n-1; cutoff <- c(0)
    if(length(Epsilon)>1) stop("Asymmetric intervals are not implemented")
    ncp <- Epsilon^2 * n; 
    findquant <- function(alpha, q, df1, df2, ncp)
        alpha - pf(q, df1, df2, ncp)
    cutoff <- sqrt(uniroot(findquant, c(0,100000000000000000000000),
                           alpha=alpha, df1=df1,
                            df2=min(df2, 4000), ncp=ncp)$root) 
    tstat <- mean/(std/sqrt(n))
    if (abs(tstat) < cutoff) result <- "rejected"
    else result <- "not rejected"
    power=2*(pt(cutoff, df2))-1
    list(Dissimilarity=result, Mean=mean, StdDev=std, n=n, alpha=alpha, 
         missing=sum(is.na(x)), Epsilon=Epsilon,
         cutoff=cutoff, Tstat=tstat, Power=power)
}

