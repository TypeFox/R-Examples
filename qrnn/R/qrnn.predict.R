qrnn.predict <-
function(x, parms)
{
    if (!is.matrix(x)) stop("\"x\" must be a matrix")
    weights <- parms$weights
    lower <- parms$lower
    eps <- min(parms$eps.seq)
    Th <- parms$Th
    x.center <- parms$x.center
    x.scale <- parms$x.scale
    y.center <- parms$y.center
    y.scale <- parms$y.scale
    lower <- (lower-y.center)/y.scale
    x <- sweep(x, 2, x.center, "-")
    x <- sweep(x, 2, x.scale, "/")
    y.bag <- matrix(0, ncol=length(weights), nrow=nrow(x))
    for (i in seq_along(weights)){
        y.bag[,i] <- qrnn.eval(x, weights[[i]]$W1, weights[[i]]$W2,
                               lower, eps, Th)
        y.bag[,i] <- y.bag[,i]*y.scale + y.center
    }
    y.bag
}
