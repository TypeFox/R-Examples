# methods.R: plotmo method functions for miscellaneous objects

plotmo.x.mars <- function(object, trace) # mda package
{
    # like plotmo.x.default but ignore object$x
    get.x.or.y(object, "x", trace, try.object.x.or.y=FALSE)
}
plotmo.type.bruto <- function(object, ...) "fitted"

plotmo.predict.bruto <- function(object, newdata, type, ..., TRACE) # mda package
{
    # TODO fails: predict.bruto returned a response of the wrong length
    plotmo.predict.defaultm(object, newdata, type=type, ..., TRACE=TRACE)
}
plotmo.type.lars <- function(object, ...) "fit"

plotmo.predict.lars <- function(object, newdata, type, ..., TRACE) # lars package
{
    # newx for predict.lars must be a matrix not a dataframe,
    # so here we use plotmo.predict.defaultm (not plotmo.predict.default)
    plotmo.predict.defaultm(object, newdata, type=type, ..., TRACE=TRACE)$fit
}
plotmo.predict.mvr <- function(object, newdata, type, ..., TRACE) # pls package
{
    # the following calls predict.mvr
    y <- plotmo.predict.default(object, newdata, type=type, ..., TRACE=TRACE)
    dim <- dim(y)
    if(length(dim) == 3) { # type="response" returns a 3 dimensional array
        if(dim[2] != 1)
            stop0("multiple response models are not supported")
        y <- y[,1,]
    }
    y
}
plotmo.predict.quantregForest <- function(object, newdata, type, ..., TRACE)
{
    # the following calls predict.quantregForest
    plotmo.predict.default(object, newdata, def.quantiles=.5, ..., TRACE=TRACE)
}
# plotmo.type.cosso works only if before calling plotmo
# you manually do class(cosso.object) <- "cosso"
plotmo.type.cosso <- function(object, ...) "fit" # cosso package

plotmo.predict.cosso <- function(object, newdata, type, ..., TRACE)
{
    # xnew for predict.cosso must be a matrix not a dataframe,
    # so here we use plotmo.predict.defaultm (not plotmo.predict.default).
    # We default M so first time users can call plotmo easily.
    yhat <- plotmo.predict.defaultm(object, newdata, type=type,
                         def.M=min(ncol(newdata), 2), ..., TRACE=TRACE)
    stopifnot(NCOL(yhat) == 1)
    # class(yhat) is "predict.cosso" but that chokes as.data.frame later
    class(yhat) <- "vector"
    yhat
}
plotmo.type.lda <- function(object, ...) "class"

plotmo.type.qda <- function(object, ...) "class"

plotmo.predict.lda <- function(object, newdata, type, ..., TRACE) # MASS package
{
    # the following calls predict.lda
    yhat <- plotmo.predict.default(object, newdata, ..., TRACE=TRACE)
    get.lda.yhat(object, yhat, type, trace=0)
}
plotmo.predict.qda <- function(object, newdata, type, ..., TRACE) # MASS package
{
    # the following calls predict.qda
    yhat <- plotmo.predict.default(object, newdata, ..., TRACE=TRACE)
    get.lda.yhat(object, yhat, type, trace=0)
}
# Special handling for MASS lda and qda predicted response, which
# is a data.frame with fields "class", "posterior", and "x".
# Here we use plotmo's type argument to choose a field.

get.lda.yhat <- function(object, yhat, type, trace)
{
    yhat1 <- switch(match.choices(type,
                        c("class", "posterior", "response", "ld"), "type"),
           class     = yhat$class, # default
           posterior = yhat$posterior,
           response  = yhat$x,
           ld        = {
                    warning0("type=\"ld\" is deprecated for lda and qda models");
                    yhat$x
                })
    if(is.null(yhat1)) {
        msg <- paste0(
            if(!is.null(yhat$x)) "type=\"response\" " else "",
            if(!is.null(yhat$class)) "type=\"class\" " else "",
            if(!is.null(yhat$posterior)) "type=\"posterior\" " else "")
        stop0("type=\"", type, "\" is not allowed for predict.",
              class(object)[1], ".  ",
              if(nzchar(msg)) paste("Use one of:", msg) else "",
              "\n")
    }
    yhat1
}
plotmo.type.varmod <- function(object, ...) "se"

plotmo.x.varmod <- function(object, trace)
{
    attr(object$parent, ".Environment") <-
        get.model.env(object$parent, "object$parent", trace)
    plotmo.x(object$parent, trace)
}
plotmo.y.varmod <- function(object, trace, naked, expected.len)
{
    attr(object$residmod, ".Environment") <-
        get.model.env(object$residmod, "object$residmod", trace)
    plotmo.y(object$residmod, trace, naked, expected.len)
}
order.randomForest.vars.on.importance <- function(object, x)
{
    importance <- object$importance
    if(!is.matrix(importance) ||    # sanity checks
       nrow(importance) == 0  ||
       !identical(row.names(importance), colnames(x))) {
        warning0("randomForest object has an invalid 'importance' field")
        return(NULL)
    }
    # vector of var indices, most important vars first
    order(importance[,1], decreasing=TRUE)
}
plotmo.singles.randomForest <- function(object, x, nresponse, trace, all1)
{
    importance <- order.randomForest.vars.on.importance(object, x)
    if(all1)
        return(importance)
    if(is.null(importance))
        seq_len(NCOL(x)) # all variables
    # 10 most important variables
    # (10 becauses plotmo.pairs returns 6, total is 16, therefore 4x4 grid)
    importance[seq_len(min(10, length(importance)))]
}
plotmo.pairs.randomForest <- function(object, x, ...)
{
    if(is.null(object$forest))
        stop0("object has no 'forest' component ",
              "(use keep.forest=TRUE in the call to randomForest)")
    importance <- order.randomForest.vars.on.importance(object, x)
    if(is.null(importance))
        return(NULL)
    # pairs of four most important variables (i.e. 6 plots)
    form.pairs(importance[1: min(4, length(importance))])
}
possible.biglm.warning <- function(object, trace)
{
    if(inherits(object, "biglm")) {
        n <- check.integer.scalar(object$n, min=1)
        y <- plotmo.y.default(object, trace, naked=TRUE, expected.len=NULL)$field
        if(NROW(y) != n)
            warnf("plotting %g cases but the model was built with %g cases\n",
                  NROW(y), n)
    }
}
plotmo.predict.biglm <- function(object, newdata, type, ..., TRACE) # biglm package
{
    # predict.biglm: newdata must include the response even though it isn't needed
    # The following extracts the response from the formula, converts it to a
    # string, then "nakens" it (converts e.g. "log(Volume)" to plain "Volume").
    resp.name <- naken(format(formula(object)[[2]]))
    if(TRACE >= 1)
        printf("plotmo.predict.biglm: adding dummy response \"%s\" to newdata\n",
               resp.name)
    data <- data.frame(NONESUCH.RESPONSE=1, newdata)
    colnames(data) <- c(resp.name, colnames(newdata))
    plotmo.predict.default(object, data, type=type, ..., TRACE=TRACE)
}
plotmo.predict.boosting <- function(object, newdata,  # adabag package
    type="prob", newmfinal=length(object$trees), ...)
{
    stopifnot(inherits(object, "boosting") || inherits(object, "bagging"))
    predict <- predict(object, newdata=newdata, newmfinal=newmfinal, ...)
    # adabag (version 4.0) returns a list, so use the type arg to select what we want
    # note that data.frames are lists, hence must check both below
    if(is.list(predict) && !is.data.frame(predict))
        predict <-
            switch(match.arg(type, c("response", "votes", "prob", "class")),
                response = predict$prob, # plotmo default, same as prob
                votes    = predict$votes,
                prob     = predict$prob,
                class    = predict$class)

    stopifnot(!is.null(predict), NROW(predict) == NROW(newdata))
    predict
}
plotmo.predict.bagging <- function(object, newdata,  # adabag package
    type="prob", newmfinal=length(object$trees), ...)
{
    plotmo.predict.boosting(object, newdata=newdata,
                            type=type, newmfinal=newmfinal, ...)
}
# TODO Following commented out because polyreg is not supported by plotmo
# So with this commented out we support plotmo(fda.object)
# but not plotmo(fda.object$fit).
# If it were not commented out, we would support neither.
#
# plotmo.singles.fda <- function(object, x, nresponse, trace, all1)
# {
#     trace2(trace, "Invoking plotmo_x for embedded fda object\n")
#     x <- plotmo_x(object$fit, trace)
#     plotmo.singles(object$fit, x, nresponse, trace, all1)
# }
# plotmo.pairs.fda <- function(object, x, nresponse, trace, all2)
# {
#     trace2(trace, "Invoking plotmo_x for embedded fda object\n")
#     x <- plotmo_x(object$fit, trace)
#     plotmo.pairs(object$fit, x, nresponse, trace, all2)
# }

# TODO following code sometimes gives bogus pairs
# plotmo.pairs.train <- function(object, x, nresponse, trace, all2) # caret package
# {
#     submod <- object$finalModel
#     # check that submod looks like an S3 model
#     if(!is.null(submod) && is.list(object))
#         plotmo.pairs(submod, x, nresponse, trace, all2)
#     else {
#         warning0("unrecognized \"train\" object") # should never get here
#         plotmo.pairs(object, x, nresponse, trace, all2)
#     }
# }

# # Simple interface for the AMORE package.
# # Thanks to Bernard Nolan and David Lorenz for these.
# # Commented out so we don't have to include AMORE in plotmo's DESCRIPTION file.
#
# plotmo.x.MLPnet <- function(object, ...)
# {
#     get("P", pos=1)
# }
# plotmo.y.MLPnet <- function(object, ...)
# {
#     get("T", pos=1)
# }
# plotmo.predict.MLPnet <- function(object, newdata, type, ..., TRACE)
# {
#     # the following calls AMORE::sim.MLPnet
#     plotmo.predict.default(object, newdata, func=AMORE::sim.MLPnet, ..., TRACE=TRACE)
# }
