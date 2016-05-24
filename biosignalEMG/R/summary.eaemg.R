summary.eaemg <- function(object, ...) {
    res <- list(empirical = object$empirical, level = object$level, gp = dim(object$intervals)[1])
    class(res) <- "summary.eaemg"
    return(res)
} 
