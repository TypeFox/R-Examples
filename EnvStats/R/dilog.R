dilog <-
function (x, ...) 
{
    if (!is.vector(x, mode = "numeric") || length(x) != 1 || 
        !is.finite(x) || x <= 1) 
        stop("'x' must be a finite numeric scalar larger than 1")
    integrate(function(y) log(y)/(y - 1), lower = 1, upper = x, 
        ...)$value
}
