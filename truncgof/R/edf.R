"edf" <-
function(x, distn = NA, parm = NA, H = NA)
{
    if (!is.na(distn) && (!is.function(try(get(distn), silent = TRUE))))
       stop("'distn' must be a character of a distribution function")
    FH <- 0
    if (is.na(H)) H <- -Inf
    if (!is.na(distn))        
       FH <- do.call(distn, c(list(H), parm))
    x <- sort(x)
    n <- length(x)
    if (n < 1) 
        stop("'x' must have 1 or more non-missing values")
    vals <- sort(unique(x))
    rval <- approxfun(vals, FH + (1-FH)*cumsum(tabulate(match(x, vals)))/n, 
        method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    attr(rval, "call") <- sys.call()
    return(rval)
}
