rdisclap <- 
function(n,p){
    if(n<0) stop("samplesize must be non negative")
    n <- floor(n)
    if (any(p < 0) | any(p >= 1))
        stop("0 <= p < 1 is required")
    if (!length(p)%in%c(1,n))
        stop("length(p) != 1 and length(p) != n")
    res <- ceiling(log(runif(n)*(1+p)/2)/log(p)-1)
    res[p == 0] <- 0
    ans <- sample(c(-1,1),n,replace=TRUE)*res
    ans <- as.integer(ans)
    return(ans)
}

