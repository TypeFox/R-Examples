`cusp.nlogLike` <-
function (p, y, X.alpha, X.beta = X.alpha, ..., verbose = FALSE) 
{
    X.alpha = as.matrix(X.alpha)
    X.beta  = as.matrix(X.beta)
    if (length(p) != NCOL(X.alpha) + NCOL(X.beta) + 2) 
        stop("p should have the as many elements as the number of columns in X.alpha plus the number\nof columns in X.beta plus 2")
    if (NCOL(y) != 1 || length(y) != NROW(X.alpha) || length(y) !=NROW(X.beta)) 
        stop("y should be vector with length equal to number of rows in X.alpha and X.beta")
    if(!exists("cusp.nc.vec")) {cusp.nc.vec <- function(...){}}
    a = p[1:NCOL(X.alpha)]
    b = p[1:NCOL(X.beta) + NCOL(X.alpha)]
    w = p[1:NCOL(y) + NCOL(X.alpha) + NCOL(X.beta)]
#    l = p[length(p) - 1]
#    s = p[length(p)]
    alpha = X.alpha %*% a
    beta = X.beta %*% b
    ab = cbind(alpha, beta)
    uab = ab[!duplicated.default((round(ab * 1000)/1000) %*% 
        c(1, 1e+05)), , drop = FALSE]
    br = sort(uab %*% c(1, 1e+05))
    mdbr <- if (length(br) > 1) 
        min(diff(br))
    else Inf
    br = c(br - mdbr/2, max(br) + mdbr)
    iab <- findInterval(ab %*% c(1, 1e+05), br)
    fab <- tabulate(iab)
    iab <- findInterval(uab %*% c(1, 1e+05), br)
    uab <- uab[order(iab), , drop = FALSE]
    if (verbose) 
#        print(cbind(alpha = uab[, 1], beta = uab[, 2], loc = l, 
#            scale = s))
		print(c(a=a, b=b, scaled.projection.coefs=w))
    z = y %*% w #(y - l)/abs(s)
    s = 1/w[NCOL(y)]
    -1 * sum(alpha * z + beta * z^2/2 - z^4/4) + 1 * sum(fab * 
        log(abs(s) * cusp.nc.vec(uab[, 1], uab[, 2])))
}

