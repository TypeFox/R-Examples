# Vaguely adapted from bootstrap package by Rob Tibshirani (R port by Friedrich Leisch)
# License: BSD
#
# Three main modifications:
# 1) Estimate the model only once for each cell of the table,
# and compute a mean weighted by cell frequencies at the end
# (Wong, Association models, 2010, p. 28-29; Clogg & Shihadeh, 1994, p. 34-38)
# 2) theta must return a vector rather than a single number,
# to allow handling a series of parameters in the same run
# 3) The variance-covariance matrix is returned, instead of standard errors
#
# Formula taken from original function and adapted to work with table cell weights.
#
# Objects needed by theta() on the cluster nodes need to be passed as arguments to
# theta().
jackknife <- function(x, theta, ...,
                      w=rep(1, length(x)),
                      cl=NULL, load.balancing=FALSE)
{
    call <- match.call()
    stopifnot(length(w) == length(x))
    w <- as.numeric(w)
    stopifnot(all(w >= 0))
    n <- length(x)
    u <- vector("list", n)

    # Run this first to find out caller errors before running parLapply
    thetahat <- as.numeric(theta(x, ...))

    if(!is.null(cl) && requireNamespace("parallel")) {
        if(load.balancing)
            u <- parallel::parLapplyLB(cl, 1:n, function(i, x, theta, ...) theta(x[-i], ...), x, theta, ...)
        else
            u <- parallel::parLapply(cl, 1:n, function(i, x, theta, ...) theta(x[-i], ...), x, theta, ...)
    }
    else {
        u <- lapply(1:n, function(i) as.numeric(theta(x[-i], ...)))
    }

    u <- do.call(rbind, u)

    # Remove replicates with NAs when computing statistics
    # This can be used by theta() to skip a failed replicate
    u2 <- u[rowSums(is.na(u)) == 0,]

    tot <- sum(w)
    mean.u <- colSums(sweep(u2, 1, w, "*"))/tot
    jack.bias <- (tot - 1) * (mean.u - thetahat)
    dev.u <- sweep(u2, 2, mean.u, "-")

    list(dev = dev.u, bias = jack.bias, values = u)
}

