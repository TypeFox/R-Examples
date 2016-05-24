rpjack <- function(x, ...){
    stopifnot(is.numeric(x))

    if(is.array(x)){
        return(apply(x, 2L:length(dim(x)), rpjack, ...))
    } else {
        x   <- as.vector(x)
        x   <- x[is.finite(x)]
        num <- length(x)
        X <- matrix(data = rep(x, num)[-(num*(0L:(num - 1L)) + seq_len(num))],
                    ncol = num,
                    nrow = num - 1L)
        return(sd(ruinprob(x = X, ...)) * (num - 1.0) / sqrt(num))
    }
}
