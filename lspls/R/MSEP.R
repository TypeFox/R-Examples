### MSEP.R: MSEP, RMSEP and R^2 functions.
### $Id: MSEP.R 40 2009-07-20 16:37:35Z bhm $


## MSEP takes a CV-object, and calculates the MSEP
MSEP.lsplsCv <- function(object, scale = FALSE, ...) {
    if (is.null(object$mode))
        stop("`object' has no `model' component.  Recalculate with `model = TRUE'")
    resp <- as.matrix(model.response(model.frame(object)))
    pred <- object$pred
    if (isTRUE(scale)) {
        sds <- sd(resp)
        resp <- sweep(resp, 2, sds, "/")
        pred <- sweep(pred, 2, sds, "/")
    }
    colMeans((pred - c(resp))^2)
}


## RMSEP is a wrapper around MSEP that returns its square root.
RMSEP.lsplsCv <- function(object, scale = FALSE, ...)
    sqrt(MSEP(object, scale, ...))


## R2 takes a CV-ojbect, and calculates the R^2
R2.lsplsCv <- function(object, ...) {
    if (is.null(object$mode))
        stop("`object' has no `model' component.  Recalculate with `model = TRUE'")
    resp <- as.matrix(model.response(model.frame(object)))
    pred <- object$pred
    SST <- apply(resp, 2, var) * (nrow(resp) - 1)
    1 - colSums((pred - c(resp))^2) / SST
}
