is.frbe <- function(x) {
    return(inherits(x, 'frbe') && 
           is.list(x) &&
           is.data.frame(x$forecasts) &&
           is.data.frame(x$features) &&
           is.vector(x$weights) &&
           is.vector(x$mean))
}
