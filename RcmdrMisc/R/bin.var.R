# bin a numeric variable

# Author: Dan Putler (revision by J. Fox, 5 Dec 04 & 5 Mar 13)

# last modified 2014-08-04

bin.var <- function (x, bins=4, method=c("intervals", "proportions", "natural"), labels=FALSE){
    method <- match.arg(method)
    
    if(length(x) < bins) {
        stop("The number of bins exceeds the number of data values")
    }
    x <- if(method == "intervals") cut(x, bins, labels=labels)
    else if (method == "proportions") cut(x, quantile(x, probs=seq(0,1,1/bins), na.rm=TRUE),
        include.lowest = TRUE, labels=labels)
    else {
        xx <- na.omit(x)
        breaks <- c(-Inf, tapply(xx, KMeans(xx, bins)$cluster, max))
        cut(x, breaks, labels=labels)
    }
    as.factor(x)
}
