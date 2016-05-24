simulate.mmglmlong1 <- function (object, nsim=1, seed=NULL, ...){
    if (nsim!=1)
        warning("Argument 'nsim' is redundant, determined by nrow(Xdesign)")
    if (!is.null(seed)) set.seed(seed)
    tmp <- table(object$longitude)
    if (min(tmp)!=max(tmp))
        stop("All subjects must have the same number of obervations.")
    N <- length(tmp)
    #   note that subnms will be character
    subnms <- names(tmp)
    n <- nrow(object$Xdesign)/N
    subobject <- mmglm1(NULL, object$Pi, object$delta, object$glmfamily,
                        object$beta, NULL, sigma=object$sigma,
                        nonstat=object$nonstat, size=object$size, msg=FALSE)
    y <- rep(NA, N*n)
    for (i in 1:N){
        subject <- subnms[i]
        subobject$Xdesign <- object$Xdesign[(object$longitude==subject),]
        #    inefficient to concatenate
        #  y <- c(y, simulate.mmglm1(subobject)$y)
        y[((i-1)*n+1):(i*n)] <- simulate.mmglm1(subobject)$y
    }
    object$y <- y
    return(object)
}

