# pint.R: plotmo functions for confidence and prediction intervals

# Handle plotmo's "level" argument.  Return a prediction interval dataframe
# with either or both of the following sets of columns.  What columns get
# returned depends on the capabilities of the object's predict method.
# For example, predict.lm allows us to return both i and ii,  and for
# earth models we can return only i.
#
#  (i)  lwr, upr               intervals for prediction of new data
#
#  (ii) cint.lwr, cint.upr     intervals for prediction of mean response

plotmo_pint <- function(object, newdata, type, level, trace, ipred, inverse.func)
{
    if(!is.specified(level))
        return(NULL)
    trace2(trace, "plotmo_pint for \"%s\" object\n", class(object)[1])
    stopifnot.string(type)
    # call plotmo.pint.xxx where xxx is object's class
    intervals <- plotmo.pint(object, newdata, type, level, trace)
    if(!is.null(intervals$lwr)) {
        intervals$lwr <-
            apply.inverse.func(inverse.func, intervals$lwr, object, trace)
        intervals$upr <-
            apply.inverse.func(inverse.func, intervals$upr, object, trace)
    }
    if(!is.null(intervals$cint.lwr)) {
        intervals$cint.lwr <-
            apply.inverse.func(inverse.func, intervals$cint.lwr, object, trace)
        intervals$cint.upr <-
            apply.inverse.func(inverse.func, intervals$cint.upr, object, trace)
    }
    print_summary(intervals, "prediction intervals", trace)
    intervals
}
# Return a data.frame with either or both of the following variables:
#  (i)  lwr, upr               intervals for prediction of new data
#  (ii) cint.lwr, cint.upr     intervals for prediction of mean response

plotmo.pint <- function(object, newdata, type, level, trace)
{
    UseMethod("plotmo.pint")
}
plotmo.pint.default <- function(object, newdata, type, level, trace)
{
    stop0("the level argument is not supported for ",
          class(object)[1], " objects")
}
plotmo.pint.lm <- function(object, newdata, type, level, trace)
{
    # lm objects with weights do not support confidence intervals on new data
    if(!is.null(object$weights))
        stop0("the level argument is not supported on lm objects with weights")
    pints <- predict(object, newdata, interval="prediction", level=level)
    cints <- predict(object, newdata, interval="confidence", level=level)
    data.frame(
        lwr      = pints[,"lwr"], # intervals for prediction of new data
        upr      = pints[,"upr"],
        cint.lwr = cints[,"lwr"], # intervals for prediction of mean response
        cint.upr = cints[,"upr"])
}
plotmo.pint.glm <- function(object, newdata, type, level, ...)
{
    if(!is.null(object$weights) && !all(object$weights == object$weights[1]))
        warnf(
"the level argument may not work correctly on glm objects built with weights")

    quant <- 1 - (1 - level) / 2 # .95 becomes .975

    predict <- predict(object, newdata, type=type, se.fit=TRUE)

    # special handling for where user used gam::gam instead of mgcv::gam
    if(class(predict)[1] == "numeric" &&
            "package:gam" %in% search()) {
        cat("\n")
        stop0("gam objects in the 'gam' package do not support ",
              "confidence intervals on new data")
    }
    data.frame(cint.lwr = predict$fit - quant * predict$se.fit,
               cint.upr = predict$fit + quant * predict$se.fit)
}
plotmo.pint.gam <- function(object, newdata, type, level, ...)
{
    plotmo.pint.glm(object, newdata, type, level)
}
plotmo.pint.quantregForest <- function(object, newdata, type, level, ...)
{
    q0 <- (1 - level) / 2   # .95 becomes .025
    q1 <- 1 - q0            # .975

    predict <- predict(object, newdata, quantiles=c(q0, q1))

    data.frame(lwr = predict[,1], upr = predict[,2])
}
plotmo.pint.earth <- function(object, newdata, type, level, ...)
{
    pints <- predict(object, newdata=newdata, type=type, interval="pint", level=level)
    if(is.null(newdata)) {
        cints <- predict(object, newdata=NULL, type=type, interval="cint", level=level)
        pints$cint.upr <- cints$upr
        pints$cint.lwr <- cints$lwr
    }
    pints
}
