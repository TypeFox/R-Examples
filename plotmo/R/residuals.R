# residuals.R: plotmo functions for residuals (the residuals, their scale, and name)

# "rinfo" is "residual info"
plotmo_rinfo <- function(object, type=NULL, residtype=type, nresponse=1,
                         standardize=FALSE, delever=FALSE, trace=0,
                         leverage.msg="returned as NA",
                         expected.levs=NULL,
                         labels.id=NULL, ...)
{
    trace2(trace,
        "----plotmo_rinfo: plotmo_resids(object, type=\"%s\", nresponse=%s)\n",
        type,
        if(is.na(nresponse)) "NA" else if(is.null(nresponse)) "NULL" else paste(nresponse))

    # TODO e.g. earth pclass nresp=1, plotmo_y returns pclass1st 0 or 1 but predict is 1, 2, 3
    if(!is.na(pmatch(type, "class")))
        stop0("plotres does not support type=\"class\" for \"%s\" objects",
              class(object)[1])

    # try calling residuals() directly
    tracex <- if(trace == 1) 0 else trace # already printed call to residuals in plotmo_meta
    plotmo_resids <- plotmo_resids(object, type, residtype,
                                   nresponse=nresponse, trace=tracex, ...)
    if(!is.null(plotmo_resids)) {
        resids <- plotmo_resids$resids
        labs   <- plotmo_resids$labs
        fitted <- plotmo_fitted(object, trace, nresponse, type, ...)$fitted
    } else {
        # trace=2 not 1 because we have already printed this message info in plotmo_meta
        if(trace >= 2)
                printf("calling predict because residuals was unsuccessful\n")
        fitted <- plotmo_predict(object, newdata=NULL, nresponse,
                    type, expected.levs, trace, inverse.func=NULL, ...)$yhat
        labs <- rownames(fitted)
        check.numeric.scalar(nresponse) # nresponse should be specified by now
        if(nresponse == 1)
            plotmo_y <- plotmo_y(object, nresponse, trace, nrow(fitted), object$levels)
        else {
            # TODO needed for e.g. rpart and lars where y has one col but predict has multiple cols
            tracex <- if(trace <= 0) -1 else trace # prevent msg in plotmo_nresponse, see note there
            plotmo_y <- try(plotmo_y(object, nresponse, tracex, nrow(fitted), object$levels),
                            silent=trace == 0)
            if(is.try.err(plotmo_y)) {
                if(trace >= 1)
                    printf(
"the call to plotmo_y was unsuccessful with nresponse=%g, trying again with nresponse=1\n",
                        nresponse)
                nresponse <- 1
                plotmo_y <- plotmo_y(object, nresponse, trace, nrow(fitted), object$levels)
                if(trace >= 1)
                    printf("plotmo_y is ok with nresponse forced to 1\n")
            }
        }
        y <- plotmo_y$y
        resids <- y - fitted
        colnames(resids) <- "resids"
        # TODO following will sometimes give the wrong results?
        if(!is.null(nresponse) && nresponse > NCOL(resids)) {
            if(trace >= 1)
                printf(
"forcing nresponse %g to 1 because response - fitted has one column\n", nresponse)
            nresponse <- 1
        }
        resids <- process.y(resids, object, type, nresponse,
                            expected.len=nrow(fitted),
                            expected.levs=expected.levs, trace, "residuals")$y
        trace2(trace,
            "generated the residuals using plotmo_predict() and plotmo_y()\n")
    }
    scale <- get.resid.scale(object, resids,
                             standardize, delever, trace, leverage.msg)
    trace2(trace, "----plotmo_rinfo: done\n")

    if(!is.null(labels.id)) # user specified labels.id?
        labs <- repl(paste(labels.id), length(resids)) # recycle if necessary

    list(resids = resids,      # numeric vector, standardize and delever not applied
         labs   = labs,        # resids names, may be NULL
         fitted = fitted,      # predicted values for newdata=NULL and given type
         scale  = scale$scale, # vector of 1s unless standardize or delever set
         name   = scale$name)  # "Residual" or "Delevered Residual" etc.
}
# return NULL if call to residuals failed
plotmo_resids <- function(object, type, residtype, nresponse, trace, ...)
{
    stopifnot.string(type)
    stopifnot.string(residtype)
    resids <- try(call.dots(stats::residuals, DROP="*", KEEP="PREFIX",
                            # following prevents reprint of residuals msg if fail
                            TRACE=if(trace == 0) -1 else trace,
                            force.object=object, force.type=residtype, ...),
                  silent=trace <= 1)
    # is.null check is for residuals(glmnet) which silently returns NULL
    if(is.try.err(resids) || is.null(resids))
        return(NULL)
    if(trace >= 2)
        print_summary(resids, "residuals is ",
                      details=if(trace>=2) 2 else -1)
    list(resids = process.y(resids, object, type, nresponse,
                    expected.len=NULL, expected.levs=NULL,
                    trace, "residuals")$y,
        labs=if(!is.null(names(resids))) names(resids) else rownames(resids))
}
get.resid.scale <- function(object, resids,
                            standardize, delever, trace, leverage.msg)
{
    scale <- repl(1, length(resids))
    name <- "Residual"
    standardize <- check.boolean(standardize)
    if(standardize) {
        scale <- plotmo_standardizescale(object)
        name <- "Standardized Residual"
    }
    delever <- check.boolean(delever)
    if(delever) {
        if(standardize) # don't allow double denormalization
            stop0("the standardize and delever arguments cannot both be set")
        hatvalues <- hatvalues1(object, "'delever'")
        hat1 <- which(hatvalues == 1)
        if(trace >= 0 && length(hat1) > 0)
            warnf("response[%s] has a leverage of one and will be %s",
                  paste.c(hat1), leverage.msg)
        scale <- 1 / sqrt(1 - hatvalues)
        name <- "Delevered Residual"
    }
    # leverages of 1 cause an inf scale, change to NA for easier handling later
    scale[is.infinite(scale)] <- NA
    check.vec(scale, "scale", length(resids), na.ok=TRUE)
    check(scale, "scale", "non-positive value", function(x) { x <= 0 }, na.ok=TRUE)
    list(scale = scale,
         name  = name)
}
# scale for standardization, inf if leverage is 1
plotmo_standardizescale <- function(object)
{
    if(inherits(object, "earth")) {
        if(is.null(object$varmod))
            stop0("\"standardize\" is not allowed because\n",
                  "the model was not built with varmod.method")
        se <- predict(object, type="earth", interval="se")
    } else if(inherits(object, "rlm"))
        se <- object$s
    else if(inherits(object, "glm"))
        se <- sqrt(summary(object)$dispersion)
    else if(inherits(object, "lm"))
        se <- sqrt(deviance(object) / df.residual(object))
    else
        stop0("'standardize' is not yet supported for this object")

    stopifnot(is.numeric(se))
    stopifnot(all(!is.na(se)), all(se > 0))

    1 / (se * sqrt(1 - hatvalues1(object, "'standardize'")))
}
hatvalues1 <- function(object, argname) # try hatvalues, specific err msg if fails
{
    hatvalues <- try(hatvalues(object))
    if(is.try.err(hatvalues))
        stop0(argname, " is not supported for this object ",
             "(the call to hatvalues failed)")
    hatvalues
}
