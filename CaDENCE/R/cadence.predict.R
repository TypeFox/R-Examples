cadence.predict <-
function(x, fit)
{
    if(!is.matrix(x)) stop("\"x\" must be a matrix")
    if("W1" %in% names(fit)) fit <- list(fit=fit)
    pred <- list()
    for(i in seq_along(fit)){
        nh <- names(fit)[i]
        hidden.fcn <- attr(fit[[nh]], "hidden.fcn")
        distribution <- attr(fit[[nh]], "distribution")
        x.center <- attr(fit[[nh]], "x.center")
        x.scale <- attr(fit[[nh]], "x.scale")
        if(attr(fit[[nh]], "stationary")) x <- x[,1,drop=FALSE]
        x.pred <- sweep(x, 2, x.center, "-")
        x.pred <- sweep(x.pred, 2, x.scale, "/")
        pred[[nh]] <- cadence.evaluate(x.pred, fit[[nh]]$W1, fit[[nh]]$W2,
                                       hidden.fcn, distribution)
    }
    if(i==1) pred <- pred[[1]]
    pred
}

