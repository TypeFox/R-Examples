rexp_mark <- function(ti, data, params){
    x <- rexp(n=1, params)
    names(x) <- "magnitude"
    return(x)
}
