simulate.dthmm <- function (object, nsim=1, seed=NULL, ...){
    if (!is.null(seed)) set.seed(seed)
    #   is inefficient having two loops 1:nsim
    y <- simulate(mchain(NULL, object$Pi, object$delta), nsim)$mc
    #------------------------
    rname <- paste("r", object$distn, sep="")
    x <- rep(NA, nsim)
    for (i in 1:nsim) {
        x[i] <- do.call(rname, c(list(n=1), getj(object$pm, y[i]),
                        getj(object$pn, i)))
    }
    object$x <- x
    object$y <- y
    return(object)
}

