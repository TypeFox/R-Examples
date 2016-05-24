`cusp.nlogLike.c` <-
function (p, y, X.alpha, X.beta = X.alpha, ..., verbose = FALSE) 
{
    X.alpha = as.matrix(X.alpha)
    X.beta  = as.matrix(X.beta)
    if (length(p) != NCOL(X.alpha) + NCOL(X.beta) + NCOL(y)) 
        stop("p should have the same number of elements as the number of columns in X.alpha plus the number of columns in X.beta plus the number of columns in y")
    if (NCOL(y) != NCOL(y) || NROW(y) != NROW(X.alpha) || NROW(y) != NROW(X.beta)) 
        stop("y should be vector/matrix with the same number of rows as X.alpha and X.beta")
    if(!exists("cusp.nc.vec")) {cusp.nc.vec <- function(...){}}
    a = p[1:NCOL(X.alpha)]
    b = p[1:NCOL(X.beta) + NCOL(X.alpha)]
    w = p[1:NCOL(y) + NCOL(X.alpha) + NCOL(X.beta)]
    alpha = X.alpha %*% a
    beta = X.beta %*% b
    if (verbose) 
	print(c(a=a,b=b,scaled.projection.coefs=w))
	g = w #c(- w[2:NCOL(y)-1],1) / w[NCOL(y)]
	s = 1/w[NCOL(y)]
	z = y %*% g
    .v = -1 * sum(alpha * z + beta * z^2/2 - z^4/4) + 1 * sum(log(abs(s) * 
        cusp.nc.vec(alpha, beta)))
	if(!is.finite(.v)) {browser()}
	.v
}

