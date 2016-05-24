pattern.n <- function(x, n) {    
    if (!is.list(x)) stop("\"x\" must be a list of matrices.\n")
    if (length(x)!=length(n)) stop("The lengths of \"x\" and \"n\" must be the same.\n")
    
    fun <- function(x1, n1) {
        ## x2: a copy of x1
        x2 <- x1
        ## replace NA with 0
        x2[is.na(x1)] <-0
        ## replace not NA with the sample size
        x2[!is.na(x1)] <- n1
        x2}
    my.df <- mapply(fun, x, n, SIMPLIFY=FALSE)
    Reduce('+', my.df)
}
