# predict.earth.R

# predict.earth returns multiple columns for multiple response models

predict.earth <- function(
    object   = stop("no 'object' argument"),
    newdata  = NULL,
    type     = c("link", "response", "earth", "class", "terms"),
                        # "terms" always returns the earth not glm terms
                        # "terms" returns just the additive terms!
                        # and just the first response if more than one
    interval = "none",
    level    = .95,     # only used if interval != none
    thresh   = .5,      # only used if type="class"
    trace    = FALSE,
    ...)                # unused, for compatibility with generic predict
{
    print.returning.earth <- function(trace, object, msg)
    {
        if(trace >= 1)
            if(is.null(object$glm.list))
                cat("predict.earth: returning earth", msg, "\n")
            else
                cat("predict.earth: returning earth (not glm)", msg, "\n")
    }
    get.predicted.response <- function(object, newdata, type)
    {
        is.type.class <- FALSE
        ylevels <- object$levels
        if(type=="class") {
            is.type.class <- TRUE
            type <- "response" # we want predicted probabilities
            if(is.null(ylevels))
                ylevels <- c(FALSE, TRUE)
        }
        if(is.null(newdata)) {
            if(is.null(object$glm.list) || type=="earth") {
                print.returning.earth(trace, object, "fitted.values")
                fit <- object$fitted.values
            } else {    # glm predictions
                trace1(trace, "predict.earth: returning glm fitted.values\n")
                fit <- matrix(0, nrow=nrow(object$fitted.values),
                                 ncol=ncol(object$fitted.values))
                colnames(fit) <- colnames(object$fitted.values)
                for(i in seq_along(object$glm.list))
                    fit[,i] = predict.glm(object$glm.list[[i]], type=type)
            }
        } else { # user supplied newdata
            bx <- model.matrix.earth(object, newdata, env=parent.frame(),
                    trace=trace,
                    Callers.name="model.matrix.earth from predict.earth")
            if(trace >= 1) {
                print_summary(bx, "predict.earth: bx", trace=2)
                trace2(trace, "\n")
            }
            if(is.null(object$glm.list) || type=="earth") {
                print.returning.earth(trace, object, "predictions")
                fit <- bx %*% object$coefficients
            } else { # glm predictions
                if(trace >= 1)
                    cat("predict.earth: returning glm", type, "predictions\n")
                fit <- matrix(0, nrow=nrow(bx),
                                 ncol=ncol(object$fitted.values))
                colnames(fit) <- colnames(object$fitted.values)
                bx <- eval(bx[,-1, drop=FALSE], env) # -1 to drop intercept
                bx.data.frame <- as.data.frame(bx)
                for(i in seq_along(object$glm.list)) {
                    fit[,i] = predict.glm(object$glm.list[[i]],
                                          newdata=bx.data.frame, type=type)
                    check.nrows(nrow(bx), nrow(fit),
                                nrow(object$fitted.values), "predict.earth")
                }
            }
        }
        if(is.type.class)
            fit <- convert.predicted.response.to.class(fit, ylevels,
                                colnames(object$coefficients)[1], thresh)
        fit
    }
    # returns just enough for termplot to work
    get.terms <- function(object, newdata)
    {
        if(!is.null(object$glm.list))
            warning0("predict.earth: returning the earth (not glm) terms")
        bx <- model.matrix.earth(object, x=newdata, env=env,
                                 trace=trace, Callers.name="predict.earth")
        dirs <- object$dirs[object$selected.terms, , drop=FALSE]
        # retain only additive terms
        additive.terms <- get.degrees.per.term(dirs) == 1
        bx <- bx[, additive.terms, drop=FALSE]
        dirs <- dirs[additive.terms, , drop=FALSE]
        coefs <- object$coefficients[additive.terms, 1, drop=FALSE]
        additive.preds <- colSums(abs(dirs)) != 0
        dirs <- dirs[, additive.preds, drop=FALSE]
        var.names <- variable.names(object, use.names=TRUE)[additive.preds]
        termMat <- matrix(0, nrow=nrow(bx), ncol=ncol(dirs))
        colnames(termMat) <- var.names
        if(ncol(bx) >= 1)
            for(ipred in seq_len(ncol(dirs)))
                for(iterm in seq_len(ncol(bx)))
                    if(dirs[iterm, ipred])
                        termMat[, ipred] =
                            termMat[, ipred] + coefs[iterm] * bx[, iterm]
        termMat
    }
    #--- predict.earth starts here ---
    check.classname(object, substitute(object), "earth")
    warn.if.dots(...)
    trace <- as.numeric(check.numeric.scalar(trace, logical.ok=TRUE))
    env <- parent.frame() # the environment from which predict.earth was called
    fit <- switch(match.arg1(type, "type"),
        "link"     = get.predicted.response(object, newdata, "link"),
        "response" = get.predicted.response(object, newdata, "response"),
        "earth"    = get.predicted.response(object, newdata, "earth"),
        "class"    = get.predicted.response(object, newdata, "class"),
        "terms"    = get.terms(object, newdata))

    # check interval argument is legal
    interval <- match.choices(interval,
                    c("none", "pint", "cint", "se", "abs.residual"), "interval")
    if(interval == "none") {
        if(!missing(level))
            stop0("predict.earth: level=", level, " was specified but interval=\"none\"")
        return(fit)
    }
    # the interval argument was used
    if(is.null(object$varmod))
        stop0("no prediction intervals because ",
              "the earth model was not built with varmod.method")
    if(type=="class" || type == "terms")
        stop0("predict.earth: the interval argument is not allowed ",
              "with type=\"", type, "\"")
    if(!is.null(object$glm.list) && type != "earth")
        stop0("predict.earth: with earth-glm models, use type=\"earth\" ",
              "when using the interval argument")
    if(NCOL(fit) != 1)
        stop0("predict.earth: the interval argument is not supported ",
              "for multiple response models")
    predict.varmod(object$varmod, newdata=newdata, type=interval, level=level)
}
convert.predicted.response.to.class <- function(resp, ylevels, resp.name, thresh=.5)
{
    which1 <- function(row, thresh) # row is a scalar or a vector
    {
        if(length(row) > 1)
            which. <- which.max(row)
        else
            which. <- if(row > thresh) 2 else 1
        which.
    }
    if(is.null(ylevels)) # should never happen
        stop0("cannot use type=\"class\" with this model")
    resp <- ylevels[apply(resp, 1, which1, thresh)]
    if(is.character(ylevels))
        resp <- factor(resp, levels = ylevels)
    fit <- as.matrix(resp, ncol=1)
    colnames(fit) <- resp.name
    fit
}
