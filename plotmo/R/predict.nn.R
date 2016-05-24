# predict.nn.R: plotmo support for the neuralnet package
#               Note that the neuralnet package is not the V&R nnet package.
#
# The neuralnet function doesn't save the standard terms etc., so we
# have to do things in a slightly non-standard way below.
#
# The rep argument must be "mean" (return mean of predicted value over all
# reps) or "best" (return predicted value on best rep) or a column index
# (return predicted value from the given rep), or an integer vector
# (return mean of predicted value over the given reps)
#
# Some of the error tests below may be duplicated in neuralnet::compute,
# but we do them here just to be sure and to avoid obscure failures later,
# and also to detect if internal implementation of nn objects changes.
#
# TODO error handling in this function hasn't been completely tested

predict.nn <- function(object, newdata=NULL, rep="mean", trace=FALSE, ...)
{
    stop.if.dots(...)  # "..." is required for compat with the
                       # generic predict, although we don't use it

    stopifnot(is.numeric(trace) || is.logical(trace), length(trace) == 1)

    if(is.null(newdata))
        newdata <- object$covariate
    stopifnot(length(dim(newdata)) == 2)
    if(NCOL(newdata) != NCOL(object$covariate))
        stop0("newdata has ", NCOL(newdata), " columns but original data had ",
              NCOL(object$covariate), " columns")
    varnames <- object$model.list$variables
    if(!is.null(colnames(newdata)) && !is.null(varnames)) {
        stopifnot(length(colnames(newdata)) == length(varnames))
        if(any(colnames(newdata) != varnames))
            warning0("colnames(newdata) do not match the ",
                     "colnames of the original data\n",
                     "         colnames(newdata): ",
                     paste.trunc(colnames(newdata)), "\n",
                     "         colnames(orginal): ", paste.trunc(varnames))
    }
    check.df.numeric.or.logical(newdata)

    result.matrix <- object$result.matrix
    if(is.null(result.matrix)) {
        # following happens if neuralnet() gave warning "algorithm did not converge"
        stop0("predict.nn: object does not have a result.matrix (did neuralnet converge?)")
    }
    stopifnot(length(dim(result.matrix)) == 2)

    stopifnot(is.character(rep) || is.numeric(rep))
    reps <- rep
    if(is.character(rep))
        switch(match.choices(rep[1], c("best", "mean"), "rep"),
            best = { reps <- which.min(result.matrix["error",])
                     if(trace)
                         cat("predict.nn: rep = \"best\" is rep =",
                             reps, "\n")
                   },
            mean = { reps <- seq_len(NCOL(result.matrix))
                     if(trace)
                         cat("predict.nn: rep = \"mean\" will take the mean of",
                             length(reps), "reps\n")
                   })

    stopifnot(!is.null(reps))
    mean.yhat <- rep_len(0, NROW(newdata))
    for(rep in reps) {
        stopifnot(length(rep) == 1, floor(rep) == rep,
                  rep >= 1, rep <= NCOL(result.matrix))
        yhat <- neuralnet::compute(x=object, covariate=newdata, rep=rep)$net.result
        stopifnot(NROW(yhat) == NROW(newdata))
        mean.yhat <- mean.yhat + yhat
    }
    mean.yhat / length(reps)
}
# plotmo method for predict.nn
# this wrapper is used merely to pass trace.call.global to predict.nn

plotmo.predict.nn <- function(object, newdata, type, ..., TRACE, FUNC=NULL)
{
    # the following invokes predict.nn
    plotmo.predict.default(object, newdata, # type arg is unused
                           trace=trace.call.global >= 1, ..., TRACE=TRACE)
}
