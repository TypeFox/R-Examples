is.farules <- function(x) {
    return(inherits(x, 'farules') && 
           is.list(x) &&
           is.list(x$rules) &&
           is.matrix(x$statistics))
}
