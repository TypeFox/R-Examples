rsad <- function(S, frac, sad, Pois.samp=TRUE, k, zeroes=FALSE, ssize=1, ...){
    if(ssize<1)stop("ssize must be at least one")
	sad <- get (paste("r", sad, sep=""), mode = "function")
    dots <- list(...)
    com <- do.call(sad,c(list(n=S),dots))
    if(Pois.samp) sam=rpois(S*ssize,lambda=frac*com)
    else {
        if(missing(k))stop("For negative binomial sampling please provide a value for k")
        else sam <- rnbinom(S*ssize,mu=frac*com,size=k)
    }
    if(ssize>1){
        y <- data.frame(sample=rep(1:ssize,each=S), species=rep(1:S,ssize), abundance=sam)
        if(!zeroes) y <- y[y$abundance>0,]
    }
    else {
        y <- sam
        if(!zeroes) y <- y[y>0]
    }
    return(y)
}
