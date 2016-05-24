`rescale.symcoca` <- function(object, axes = seq_len(object$n.axes), ...) {
    col.nam <- colnames(object$scores$species$Y)[axes]
    lambda4 <- diag(sqrt(sqrt(object$lambda[axes])),
                    nrow = length(axes),
                    ncol = length(axes))
    Y.scale <- object$scores$species$Y[, axes, drop = FALSE] %*% lambda4
    X.scale <- object$scores$species$X[, axes, drop = FALSE] %*% lambda4
    colnames(Y.scale) <- colnames(X.scale) <- col.nam
    retval <- list(Y = Y.scale, X = X.scale)
    retval
}

