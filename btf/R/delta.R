##' generate Matrix D^{k+1}
##'
##' @param n sample size
##' @param k order of fit
##' @param x vector of inputs on domain
##' @author Edward A. Roualdes
genDelta <- function(n, k, x) {
    nk <- n-k-1
    d <- Matrix(0, nk, n)
    for (i in seq_len(nk)) {
        ik <- i+k+1
        idx <- i:ik
        tmp <- x[ik] - x[i]
        z <- x[idx]
        d[i,idx] <- tmp*sapply(seq_along(z),function(j) 1/prod(z[j]-z[-j]))
    }
    gamma(k+1)*d/n^k
}

