init.centers <- function(x, k, method = c("kmeans++", "random"))
{
    method.int <- switch ( match.arg(method), "kmeans++" = 0, "random" = 1)
    
    if (! is.matrix(x)) stop("init.centers: x is not a matrix")
    if (k < 0) stop("init.centers: k < 0")
    
    n <- ncol(x)
    m <- nrow(x)
    
    centers <- NULL
    
    if (method.int == 0)
        centers <- .Call(init_kmeanspp_r , x, k)
    else if (method.int == 1)
        centers <- .Call(init_random_r , x, k)
    
    centers  
}