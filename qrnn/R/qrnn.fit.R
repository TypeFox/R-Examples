qrnn.fit <-
function(x, y, n.hidden, tau=0.5, n.ensemble=1, iter.max=5000,
         n.trials=5, bag=FALSE, lower=-Inf,
         eps.seq=2^(-8:-32), Th=sigmoid, Th.prime=sigmoid.prime,
         penalty=0, trace=TRUE, ...)

{
    if (!is.matrix(x)) stop("\"x\" must be a matrix")
    if (!is.matrix(y)) stop("\"y\" must be a matrix")
    if (ncol(y) != 1) stop("\"y\" must be univariate")
    if ((tau > 1) | (tau < 0)) stop("invalid \"tau\"")
    if (identical(Th, linear)) n.hidden <- 1
    x <- scale(x)
    x.center <- attr(x, "scaled:center")
    x.scale <- attr(x, "scaled:scale")
    y <- scale(y)
    y.center <- attr(y, "scaled:center")
    y.scale <- attr(y, "scaled:scale")
    lower.scaled <- (lower-y.center)/y.scale
    weights <- list()
    if(trace) cat("tau=", tau, "\n", sep="")
    for (i in 1:n.ensemble){
        if(trace) cat(i, "/", n.ensemble, "\n", sep="")
        w.tmp <- NA
        class(w.tmp) <- "try-error"
        n.errors <- 0
        while (inherits(w.tmp, "try-error")) {
            w.tmp <- try(qrnn.nlm(x, y, n.hidden, tau, iter.max,
                                  n.trials, bag, lower.scaled, eps.seq,
                                  Th, Th.prime, penalty, trace, ...),
                        silent = TRUE)
            n.errors <- n.errors + 1
            if (n.errors > 5) stop("nlm optimization failed")
        }
        weights[[i]] <- w.tmp
    }
    if(trace) cat("\n")
    parms <- list(weights=weights, lower=lower, eps.seq=eps.seq,
                  tau=tau, Th=Th, x.center=x.center,
                  x.scale=x.scale, y.center=y.center, y.scale=y.scale)
    parms
}
