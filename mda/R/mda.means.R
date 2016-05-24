mda.means <-
function (object, x, y) 
{
    weights <- means <- object$weights
    nn <- names(object$weights)
    for (i in nn) {
        xx <- x[y == i, ]
        ww <- weights[[i]]
        nc <- ncol(ww)
        xm <- matrix(0, ncol(x), nc)
        for (j in seq(nc)) {
            www <- ww[, j]
            www <- www/sum(www)
            xm[, j] <- apply(xx * www, 2, sum)
        }
        means[[i]] <- xm
    }
    means
}

