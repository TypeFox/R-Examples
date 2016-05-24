`predict.timetrack` <- function(object, newdata, ...) {
    namNew <- deparse(substitute(newdata))
    ## Apply a transformation - let tran deal with arg matching
    if(!isTRUE(all.equal(transform, "none"))) {
        newdata <- tran(newdata, method = object$transform, ...)
    }
    ## merge X and passive
    dat <- join(object$X, newdata, type = "left")
    X <- dat[[1]]
    newdata <- dat[[2]]
    ## common set of species
    tmp <- colSums(X > 0) > 0
    X <- X[, tmp]
    newdata <- newdata[, tmp]

    ## fitted values for newdata
    pred <- predict(object$ordination, newdata = newdata, type = "wa",
                    scaling = object$scaling, model = "CCA",
                    rank = object$rank)
    pred2 <- predict(object$ordination, newdata = newdata, type = "wa",
                     scaling = object$scaling, model = "CA",
                     rank = object$rank)
    pred <- cbind(pred, pred2)
    ## return object
    nams <- object$labels
    nams[["passive"]] <- namNew
    ## update object with the new passive data predictions
    object$fitted.values <- pred
    object$labels <- nams
    object
}
