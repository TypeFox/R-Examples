number.suffix <-
function (x) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    if (any(is.na(x))) 
        stop("Missing values not allowed in 'x'")
    n <- length(x)
    ret.vec <- character(n)
    for (i in 1:n) ret.vec[i] <- number.suffix.scalar(x[i])
    ret.vec
}
