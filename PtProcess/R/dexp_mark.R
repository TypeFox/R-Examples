dexp_mark <- function(x, data, params){
    y <- dexp(x[,"magnitude"], params, log=TRUE)
    return(y)
}
