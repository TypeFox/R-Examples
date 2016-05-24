# $Id: mvnorm.R 275 2014-05-02 22:07:27Z mmaechler $

rmvnorm <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                   method=c("eigen", "svd", "chol"), pre0.9_9994 = FALSE)
{
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                     check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma))
        stop("mean and sigma have non-conforming size")

    method <- match.arg(method)

    R <- if(method == "eigen") {
        ev <- eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))){
            warning("sigma is numerically not positive definite")
        }
        ## ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
        ## faster for large  nrow(sigma):
        t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
    }
    else if(method == "svd"){
        s. <- svd(sigma)
        if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))){
            warning("sigma is numerically not positive definite")
        }
        t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
    }
    else if(method == "chol"){
        R <- chol(sigma, pivot = TRUE)
        R[, order(attr(R, "pivot"))]
    }

    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*%  R
    retval <- sweep(retval, 2, mean, "+")
    colnames(retval) <- names(mean)
    retval
}

dmvnorm <- function (x, mean = rep(0, p), sigma = diag(p), log = FALSE)
{
    if (is.vector(x))
	x <- matrix(x, ncol = length(x))
    p <- ncol(x)
    if(!missing(mean)) {
	if(!is.null(dim(mean))) dim(mean) <- NULL
	if (length(mean) != p)
	    stop("mean and sigma have non-conforming size")
    }
    if(!missing(sigma)) {
	if (p != ncol(sigma))
	    stop("x and sigma have non-conforming size")
	if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
			 check.attributes = FALSE))
	    stop("sigma must be a symmetric matrix")
    }

    ## <faster code contributed by Matteo Fasiolo mf364 at bath.ac.uk
    dec <- tryCatch(chol(sigma), error=function(e)e)
    if (inherits(dec, "error")) {
        ## warning("cannot compute chol(sigma)"); return(NaN)
        ## behave the same as dnorm(): return Inf or 0
        x.is.mu <- colSums(t(x) != mean) == 0
        logretval <- rep.int(-Inf, nrow(x))
        logretval[x.is.mu] <- Inf # and all other f(.) == 0
    } else {
	tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
	rss <- colSums(tmp ^ 2)
	logretval <- - sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
    }
    names(logretval) <- rownames(x)
    if(log) logretval else exp(logretval)
}
