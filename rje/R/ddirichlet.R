ddirichlet <-
function (x, alpha, log = FALSE, tol = 1e-10) 
{
    k = length(alpha)
    if (is.matrix(x)) 
        n = dim(x)[1]
    else if (is.numeric(x)) {
        dim(x) = c(1, length(x))
        n = 1
    }
    else stop("Vales must be numeric")
    if (any(alpha <= 0)) 
        stop("Parameters must be positive")
    out = numeric(n)
    if (dim(x)[2] == k - 1) 
        x = cbind(x, 1 - rowSums(x))
    else if (dim(x)[2] != k) 
        stop("Dimensions do not match")
    out[(apply(x, 1, function(y) any(y < -tol))) | abs(rowSums(x) - 
        1) > tol] = 0
    out = colSums((alpha - 1) * log(t(x))) + lgamma(sum(alpha)) - 
        sum(lgamma(alpha))
    if (!log) 
        out = exp(out)
    return(out)
}
