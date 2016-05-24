UnivTest <- function (x, type = c("truncated", "bartlett", "daniell", "QS",    
 "parzen"), testType = c("covariance","correlation"), p, b = 0, parallel = FALSE) 
{
    type <- match.arg(type)
    testType <- match.arg(testType)
    if (missing(testType)) method="covariance"  
    data.name <- deparse(substitute(x))
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
    
    MaxLag <- n - 1
    adcv <- function(k,x){
        n <- length(x)
        xA <- x[1:(n-k)]
        xB <- x[(1+k):n]
        adcv <- dcov(xA,xB)
        return(adcv)
      }
    adcf <- function(k,x){
        n <- length(x)
        xA <- x[1:(n-k)]
        xB <- x[(1+k):n]
        adcf <- dcor(xA,xB)
        return(adcf)
      }
    t = rep(0, MaxLag)
    for (k in 1:MaxLag) {
        kern <- kernelFun(type, k/p)
        if (kern != 0) {
                if(testType=="covariance"){
                 t[k] <- (n - k) * kern^2 * adcv(k,x)^2
                } else {
                 t[k] <- (n - k) * kern^2 * adcf(k,x)^2 
                }
        }
    }
    method2 = ifelse((testType=="covariance"),"Univariate test of independence based on distance covariance",
"Univariate test of independence based on distance correlation")
    stat <- sum(t)
    if (!b == 0) {
        Tnstar <- TstarBoot(x, type, testType, p, b, parallel)
        pvalue <- sum(Tnstar >= stat)/(b + 1)
    }
    p.value <- ifelse(b == 0, NA, pvalue)
    if (b == 0) {
        Tnstar <- NULL
    }
    else Tnstar <- Tnstar
    dataname <- paste(data.name, ",", " kernel type: ", type, 
        ", bandwidth=", p, ", boot replicates ", b, sep = "")
    names(stat) <- "Tn"
    e = list(method = method2, statistic = stat, p.value = p.value, replicates = Tnstar, 
        data.name = dataname)
    class(e) = "htest"
    return(e)
}
