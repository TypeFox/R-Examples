## Methods

"[.mboostLSS" <- function(x, i, return = TRUE, ...) {
    stopifnot((length(i) == 1 | length(i) == length(x)) && i > 0)
    attr(x, "subset")(i)
    if (return) return(x)
    invisible(NULL)
}

coef.mboostLSS <- function(object, which = NULL,
                           aggregate = c("sum", "cumsum", "none"),
                           parameter = names(object), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    RET <- lapply(parameter, function(i, object)
                  coef(object[[i]], which = which,
                       aggregate = aggregate, ...),
                  object = object)
    if (length(RET) == 1)
        RET <- RET[[1]]
    return(RET)
}

coef.glmboostLSS <- function(object, which = NULL,
                             aggregate = c("sum", "cumsum", "none"),
                             off2int = FALSE, parameter = names(object), ...){
    coef.mboostLSS(object, which = which, aggregate = aggregate,
                   parameter = parameter, off2int = off2int, ...)
}

risk.mboostLSS <- function(object, merge = FALSE, parameter = names(object), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)

    lo <- length(unique(mstop(object)))
    if (merge) {
        get_rsk <- function(i, object) {
            mmo <- max(mstop(object))
            rsk <- object[[i]]$risk()
            if (length(rsk) != mmo) {
                rsk <- c(rsk, rep(NA, mmo - length(rsk)))
            }
            rsk
        }
        RES <- sapply(parameter, get_rsk,
                      object = object)
        RES <- as.vector(t(RES))
        names(RES) <- rep(names(parameter), mstop(object)[1])
        ## drop unwanted NAs
        if (lo != 1)
            RES <- RES[!is.na(RES)]
        class(RES) <- object[[1]]$control$risk
        return(RES)
    }
    RES <- lapply(parameter, function(i, object)  object[[i]]$risk(),
                  object = object)
    class(RES) <- object[[1]]$control$risk
    return(RES)
}

mstop.mboostLSS <- function(object, parameter = names(object), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    RET <- sapply(parameter, function(i, object)  object[[i]]$mstop(),
                  object = object)
    names(RET) <- names(object)[parameter]
    if (length(RET) == 1)
        RET <- RET[[1]]
    return(RET)
}

mstop.oobag <- function(object, parameter = names(object), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    RET <- sapply(parameter, function(i, object)  which.min(object[[i]]),
                  object = object)
    names(RET) <- names(object)[parameter]
    if (length(RET) == 1)
        RET <- RET[[1]]
    return(RET)
}

selected.mboostLSS <- function(object, parameter = names(object), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    RET <- lapply(parameter, function(i, object)
                                 selected(object[[i]]),
                  object = object)
    names(RET) <- names(object)[parameter]
    if (length(RET) == 1)
        RET <- RET[[1]]
    return(RET)
}


plot.glmboostLSS <- function(x, main = names(x), parameter = names(x),
                             off2int = FALSE, ...){
    if (is.character(parameter))
        parameter <- extract_parameter(x, parameter)
    lapply(parameter, function(i, x, main, off2int,  ...)
                    plot(x[[i]], main = main[[i]], off2int = off2int,  ...),
           x = x, main = main, off2int = off2int, ...)
    invisible(coef(x, aggregate = "cumsum", off2int = off2int))
}


plot.gamboostLSS <- function(x, main = names(x), parameter = names(x), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(x, parameter)
    RET <- lapply(parameter, function(i, x, main, ...)
                      plot(x[[i]], main = main[[i]], ...),
                  x = x, main = main, ...)
    if (any(sapply(RET, class) == "trellis")) {
        return(RET)
    } else {
        invisible(RET)
    }
}

plot.predint <- function(x, main = "Marginal Prediction Interval(s)",
                         xlab = NULL, ylab = NULL, lty = c("solid", "dashed"),
                         lcol = c("black", "black"), log = "", ...) {

    pi <- attr(x, "pi")
    which <- attr(x, "which")
    rawdata <- attr(x, "rawdata")

    if (length(lty) != length(pi) + 1)
        lty <- c(lty, rep(tail(lty, 1), (length(pi) + 1) - length(lty)))
    if (length(lcol) != length(pi) + 1)
        lcol <- c(lcol, rep(tail(lcol, 1), (length(pi) + 1) - length(lcol)))

    if (is.null(xlab))
        xlab <- which
    if (is.null(ylab))
        ylab <- "prediction"

    plot(rawdata$x, rawdata$y, pch = 20,
         col = rgb(0.5, 0.5, 0.5, 0.5),
         xlab = xlab, ylab = ylab, main = main,
         log = log, ...)

    lines(x[, which], x$"Prediction (Median)",
          lty = lty[1], col = lcol[1], ...)

    for (i in seq_along(pi)) {
        lines(x[, which], x[, paste0(pi[i] * 100, "% PI (lower)")],
              lty = lty[i + 1], col = lcol[i + 1], ...)
        lines(x[, which], x[, paste0(pi[i] * 100, "% PI (upper)")],
              lty = lty[i + 1], col = lcol[i + 1], ...)
    }
}


PI <- predint <- function(x, which, pi = 0.9, newdata = NULL, ...) {
    qfun <- get_qfun(x)

    if (length(which) != 1 || !is.character(which))
        stop("Please specify the variable for the marginal prediction interval.")

    var <- get_data(x, which = which)
    if (ncol(var) > 1 || is.factor(var))
        stop("Prediction intervals are currently only implemented for ",
             "base-learners of one numeric variable")

    pred_vars <- lapply(x, extract, what = "variable.names")
    pred_vars <- unique(unlist(pred_vars))
    if ("(Intercept)" %in% pred_vars)
        pred_vars <- pred_vars[pred_vars != "(Intercept)"]

    if (is.null(newdata)) {
        tmp <- get_data(x, which = pred_vars)
        i <- grepl(which, names(tmp))
        if (sum(i) != 1)
            stop(sQuote("which"), " is misspecified")
        newdata <- data.frame(x1 = seq(min(tmp[, i]),
                                       max(tmp[, i]), length = 150))
        colnames(newdata) <- names(tmp)[i]
        newdata[, names(tmp)[!i]] <- lapply(tmp[, !i], mean_mod)
    } else {
        i <- grepl(which, names(newdata))
        if (sum(i) != 1)
            stop(sQuote("which"), "is misspecified")
        ## check if data is ok, else give a warning
        if (nrow(unique(newdata[, !i])) != 1)
            warning("All variables but", sQuote("which"), "should be constant")
    }

    newdata[, names(x)] <- predict(x, newdata = newdata, type = "response")
    newdata$"Prediction (Median)" <- do.call(qfun, args = c(p = 0.5,
        newdata[, names(x)]))

    for (i in seq_along(pi)) {
        newdata[, paste0(pi[i] * 100, "% PI (lower)")] <- do.call(qfun,
            args = c(p = (1 - pi[i])/2, newdata[, names(x)]))
        newdata[, paste0(pi[i] * 100, "% PI (upper)")] <- do.call(qfun,
            args = c(p = 1 - (1 - pi[i])/2, newdata[, names(x)]))
    }
    # drop predictions of parameters
    newdata <- newdata[, !names(newdata) %in% names(x)]

    class(newdata) <- c("predint", "data.frame")
    attr(newdata, "pi") <- pi
    attr(newdata, "which") <- which
    attr(newdata, "rawdata") <- data.frame(x = get_data(x, which = pred_vars)[, which],
                                           y = x[[1]]$response)
    return(newdata)
}


print.mboostLSS <- function(x, ...){
    cl <- match.call()
    cat("\n")
    cat("\t LSS Models fitted via Model-based Boosting\n")
    cat("\n")
    if (!is.null(attr(x, "call")))
        cat("Call:\n", deparse(attr(x, "call")), "\n\n", sep = "")
    cat("Number of boosting iterations (mstop): ",
        paste(names(mstop(x)), mstop(x), sep = " = ", collapse = ", "), "\n")
    nus <- sapply(x, function(xi) xi$control$nu)
    cat("Step size: ",
        paste(names(nus), nus, sep = " = ", collapse = ", "), "\n\n")
    cat("Families:\n")
    lapply(x, function(xi) show(xi$family))
    invisible(x)
}


fitted.mboostLSS <- function(object, parameter = names(object), ...){
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    myApply(parameter, function(i, mod, ...) fitted(mod[[i]], ...),
           mod = object, ...)
}


predict.mboostLSS <- function(object, newdata = NULL,
                              type = c("link", "response", "class"),
                              which = NULL,
                              aggregate = c("sum", "cumsum", "none"),
                              parameter = names(object), ...) {
    if (is.character(parameter))
        parameter <- extract_parameter(object, parameter)
    myApply(parameter, function(i, mod, ...)
            predict(mod[[i]], newdata = newdata, type = type, which = which,
                    aggregate = aggregate, ...),
            mod = object, ...)
}

update.mboostLSS <- function(object, weights, oobweights = NULL,
                             risk = NULL, trace = NULL, mstop = NULL, ...) {
    attr(object, "update")(weights = weights, oobweights = oobweights,
                           risk = risk, trace = trace, mstop = mstop, ...)
}

## generic version of model.weights (see stats::model.weights)
model.weights <- function(x, ...)
    UseMethod("model.weights")

model.weights.default <- function(x, ...)
    stats::model.weights(x)

model.weights.mboostLSS <- function(x, ...)
    attr(x, "(weights)")

### summary function based on print.mboostLSS() and summary.mboost()
summary.mboostLSS <- function(object, ...) {
    cat("\n")
    cat("\t LSS Models fitted via Model-based Boosting\n")
    cat("\n")
    if (!is.null(attr(object, "call")))
        cat("Call:\n", deparse(attr(object, "call")), "\n\n", sep = "")
    cat("Number of boosting iterations (mstop): ",
        paste(names(mstop(object)), mstop(object),
              sep = " = ", collapse = ", "),
        "\n")
    nus <- sapply(object, function(xi) xi$control$nu)
    cat("Step size: ",
        paste(names(nus), nus, sep = " = ", collapse = ", "), "\n\n")

    cat("Families:\n")
    lapply(object, function(xi) show(xi$family))

    if (inherits(object, "glmboostLSS")) {
        cat("Coefficients:\n")
        cf <- coef(object, off2int = TRUE)
        for (i in 1:length(cf)) {
            cat("Parameter ", names(cf)[i], ":\n", sep = "")
            print(cf[[i]])
            cat("\n")
        }
    }

    cat("Selection frequencies:\n")
    for (i in 1:length(object)) {
        cat("Parameter ", names(object)[i], ":\n", sep = "")
        nm <- variable.names(object[[i]])
        selprob <- tabulate(selected(object[[i]]), nbins = length(nm)) /
            length(selected(object[[i]]))
        names(selprob) <- names(nm)
        selprob <- sort(selprob, decreasing = TRUE)
        selprob <- selprob[selprob > 0]
        print(selprob)
    }
    invisible(object)
}

################################################################################
### helpers

## extract parameter index from mboostLSS object x
extract_parameter <- function(x, parameter) {
    idx <- sapply(parameter, function(w) {
        wi <- grep(w, names(x), fixed = TRUE)
        if (length(wi) > 0) return(wi)
        return(NA)
    })
    if (any(is.na(idx)))
        warning(paste(parameter[is.na(idx)], collapse = ","), " not found")
    parameter <- idx
}

## function for weighted sd
weighted.sd <- function(x, w, ...) {
    if (missing(w))
        w <- rep(1, length(x))
    m <- weighted.mean(x, w, ...)
    var <- weighted.mean((x - m)^2, w, ...) * sum(w) / (sum(w) - 1)
    return(sqrt(var))
}

## weighted median
weighted.median <- function (x, w = 1, na.rm = FALSE) {
    if (length(w) == 1)
        w <- rep(w, length(x))

    ## remove observations with zero weights
    x <- x[w != 0]
    w <- w[w != 0]

    ## remove NAs if na.rm = TRUE
    if (na.rm) {
        keep <- !is.na(x) & !is.na(w)
        x <- x[keep]
        w <- w[keep]
    } else {
        if (any(is.na(x)) | any(is.na(w)))
            return(NA)
    }

    ## sort data and weights
    ind <- order(x)
    x <- x[ind]
    w <- w[ind]

    ## first time that fraction of weights is above 0.5
    ind1 <- min(which(cumsum(w)/sum(w) > 0.5))

    ## first time that fraction of weights is below 0.5
    ind2 <- ifelse(ind1 == 1, 1, max(which(cumsum(w)/sum(w) <= 0.5)))

    ## if sum of weights is an even integer
    if(sum(w) %% 1 == 0 && sum(w) %% 2 == 0)
        return(mean(c(x[ind1], x[ind2])))

    ## else return
    return(max(c(x[ind1], x[ind2])))
}
