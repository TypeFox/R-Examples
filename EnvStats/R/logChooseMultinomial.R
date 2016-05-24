logChooseMultinomial <-
function (n, m) 
{
    if (length(n) > 1) 
        stop("n must have length 1")
    if (n != round(n)) 
        stop("n must be an integer")
    if (any(m != round(m) | m < 0)) 
        stop("Values of m must be positive integers")
    if (sum(m) != n) 
        stop("Must have sum(m) == n")
    return(lgamma(n + 1) - sum(lgamma(m + 1)))
}
