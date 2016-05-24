simulate.mchain <- function(object, nsim=1, seed=NULL, ...){
    #    simulate a Markov Chain
    if (!is.null(seed)) set.seed(seed)
    m <- ncol(object$Pi)
    #------------------------
    if (sum(object$delta)!=1) stop("Invalid delta")
    if (any(object$delta==1))
        initial <- (1:m)[as.logical(object$delta)]
    else
        initial <- sample(m, 1, prob=object$delta)
    #------------------------
    x <- rep(NA, nsim)
    x[1] <- initial
    for (i in 2:nsim)
        x[i] <- sample(x=1:m, size=1, prob=object$Pi[(x[i-1]),])
    object$mc <- x
    return(object)
}

