ADCF <- function (x, MaxLag=15) {
    n <- length(x)
    if (is.matrix(x)) {
        if (!NCOL(x) == 1) 
            stop("Univariate time series only")
    }
    else {
        x <- c(x)
    }
    if (!is.numeric(x)) 
        stop("'x' must be numeric")
    if (!all(is.finite(x))) 
        stop("Missing or infitive values")
    if((MaxLag<0) | (MaxLag>(n-1)))stop("'MaxLag' must be in the range of 0 and (n-1)")
    adcf <- matrix(0, 1,MaxLag+1,dimnames=list("  ",paste("lag",0:MaxLag)))
    for (k in seq(0,MaxLag,by=1)){
     xA <- x[1:(n-k)]
     xB <- x[(1+k):n]
     adcf[,(k+1)] <- round(dcor(xA,xB),4)
    }
    return(adcf)
}

