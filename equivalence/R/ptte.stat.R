# $Id: ptte.stat.R,v 1.2 2005/09/26 02:29:00 andrewr Exp $

"ptte.stat" <-
function(mean, std, n, alpha=0.05, Epsilon=0.25) {
    df1 <- 1; df2 <- n-1; cutoff <- c(0)
    if(length(Epsilon)>1) stop("Asymmetric intervals are not yet implemented")
    ncp <- Epsilon^2 * n; 
    findquant <- function(alpha, q, df1, df2, ncp)
        alpha - pf(q, df1, df2, ncp)
    cutoff <- sqrt(uniroot(findquant, c(0,1000000000000000000000000),
                           alpha=alpha, df1=df1,
                            df2=min(df2, 1000), ncp=ncp)$root)
    tstat <- mean/(std/sqrt(n))
    if (abs(tstat) < cutoff) result <- "rejected"
    else result <- "not rejected"
    power=2*(pt(cutoff, df2))-1
    list(Dissimilarity=result, Mean=mean, StdDev=std, n=n, alpha=alpha, 
         Epsilon=Epsilon, cutoff=cutoff, Tstat=tstat, Power=power)
}

