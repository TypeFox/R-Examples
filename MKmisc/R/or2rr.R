## odds ratio (OR) to relative risk (RR)
or2rr <- function(or, p0){
    if(any(or <= 0))
        stop("'or' has to be positive")
    if(p0 <= 0 | p0 >= 1)
        stop("'p0' has to be in (0,1)")
    if(length(p0) != 1)
        stop("'p0' has to be of length 1")
    names(or) <- NULL
    or/(1 - p0 + p0*or)
}
