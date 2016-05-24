linCombDist <- function(x, distfun1, w1, distfun2, w2, diag = FALSE, upper = FALSE){ 
    stopifnot(is.function(distfun1))
    stopifnot(is.function(distfun2))
    
    res <- w1*distfun1(x) + w2*distfun2(x)
    attributes(res) <- NULL
    attr(res, "Size") <- nrow(as.matrix(x))
    attr(res, "Labels") <- dimnames(x)[[1L]]
    attr(res, "Diag") <- diag
    attr(res, "Upper") <- upper
    attr(res, "call") <- match.call()
    class(res) <- "dist"
    res
}
