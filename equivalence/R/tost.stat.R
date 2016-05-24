# $Id: tost.stat.R,v 1.2 2005/09/26 02:29:00 andrewr Exp $

"tost.stat" <-
function(mean, std, n, null=0, alpha=0.05, Epsilon=0.36) {
    df <- n-1
    if(length(Epsilon)>1) stop("Asymmetric intervals are not yet implemented")
    tint <- (std/sqrt(n)) * qt(1-alpha, df)
    if ((mean - null + tint) < Epsilon & (mean - null - tint) > -Epsilon)
      result <- "rejected"
    else result <- "not rejected"
    out <- list(Dissimilarity=result, Mean=mean, StdDev=std, n=n, alpha=alpha, 
         Epsilon=Epsilon, Interval = as.numeric(formatC(tint, format="f",
                            digits=3)))
    return(out)
}

