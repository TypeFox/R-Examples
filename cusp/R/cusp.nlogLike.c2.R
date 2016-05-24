`cusp.nlogLike.c2` <-
function (p, y, X.alpha, X.beta = X.alpha, ..., verbose = FALSE) 
{
    X.alpha = as.matrix(X.alpha)
    X.beta  = as.matrix(X.beta)
    if (length(p) != NCOL(X.alpha) + NCOL(X.beta) + 2) 
        stop("p should have the same number of elements as the number of columns in X.alpha plus the number of columns in X.beta plus 2")
    if (NCOL(y) != 1 || length(y) != NROW(X.alpha) || length(y) != NROW(X.beta)) 
        stop("y should be vector with length equal to number of rows in X.alpha and X.beta")
    if(!exists("cusp.nc.vec")) {cusp.nc.vec <- function(...){}}
    a = p[1:NCOL(X.alpha)]
    b = p[1:NCOL(X.beta) + NCOL(X.alpha)]
    l = p[length(p) - 1]
    s = p[length(p)]
    alpha = X.alpha %*% a
    beta = X.beta %*% b
    if (verbose) 
#        print(cbind(alpha = alpha, beta = beta, loc = l, scale = s))
print(c(a=a,b=b,loc=l,scale=s))
    z = (y - l)/abs(s)
    -1 * sum(alpha * z + beta * z^2/2 - z^4/4) + 1 * sum(log(abs(s) * 
        cusp.nc.vec(alpha, beta)))
}

