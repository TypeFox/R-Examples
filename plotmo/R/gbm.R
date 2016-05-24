# gbm.R: plotmo functions for gbm objects
#
# TODO Add support for plotmo's level argument (quantile regression).

plotmo.prolog.gbm <- function(object, object.name, trace, ...) # invoked when plotmo starts
{
    if(is.null(object$data)) # TODO could do more if object had a call component
        stop0("object$data is NULL, ",
              "(use keep.data=TRUE in the call to gbm)")

    # "importance" is a vector of variable indices (column numbers in x), most
    # important vars first, no variables with relative.influence < 1%.  We attach
    # it to the object to avoid calling summary.gbm twice (it's expensive).

    attr(object, "importance") <- order.gbm.vars.on.importance(object)

    object
}
order.gbm.vars.on.importance <- function(object)
{
    # order=FALSE so importances correspond to orig variable indices
    importance <- summary(object, plotit=FALSE,     # calls summary.gbm
                          order=FALSE, normalize=TRUE)$rel.inf
    # NA assignment below so order() drops vars with importance < .01
    importance[importance < .01] <- NA
    importance <- order(importance, decreasing=TRUE, na.last=NA)
    # return a vector of variable indices, most important vars first
    importance[!is.na(importance)]
}
plotmo.singles.gbm <- function(object, x, nresponse, trace, all1)
{
    importance <- attr(object, "importance")
    stopifnot(!is.null(importance)) # uninitialized?
    if(all1)
        return(importance)
    # indices of vars with importance >= 1%, max of 10 variables
    # (10 becauses plotmo.pairs returns 6, total is 16, therefore 4x4 grid)
    importance[1: min(10, length(importance))]
}
plotmo.pairs.gbm <- function(object, ...)
{
    # pairs of four most important variables (i.e. 6 plots)
    importance <- attr(object, "importance")
    stopifnot(!is.null(importance)) # uninitialized?
    form.pairs(importance[1: min(4, length(importance))])
}
plotmo.x.gbm <- function(object, ...)
{
    # Return the first ntrain rows of the x matrix.  The x matrix is stored
    # with the gbm object as a vector, so we must convert it back to
    # a data.frame here, one column for each variable.
    # The use of as.integer here is copied from gbm.R.

    ntrain <- as.integer(object$train.fraction * nrow(object$data$x.order))
    x <- matrix(object$data$x, ncol=ncol(object$data$x.order))
    x <- data.frame(x[seq_len(ntrain), ])
    colnames(x) <- colnames(object$data$x.order)

    # convert numeric columns that are actually factors
    # TODO this only works correctly if default ordering of factors was used

    for(i in seq_len(ncol(x)))
        if(typeof(object$var.levels[[i]]) == "character")
            x[[i]] <- factor(x[[i]], labels=object$var.levels[[i]])
    x
}
plotmo.y.gbm <- function(object, ...)
{
    ntrain <- as.integer(object$train.fraction * nrow(object$data$x.order))
    object$data$y[seq_len(ntrain)]
}
plotmo.predict.gbm <- function(object, newdata, type, ..., TRACE)
{
    # The following invokes predict.gbm.
    # n.trees is defaulted so first time users can call plotmo(gbm.model) easily.
    # predict.gbm doesn't do partial matching on type so we do it here with pmatch.
    plotmo.predict.default(object, newdata,
        type = match.choices(type, c("link", "response"), "type"),
        def.n.trees = object$n.trees, ..., TRACE=TRACE)
}
