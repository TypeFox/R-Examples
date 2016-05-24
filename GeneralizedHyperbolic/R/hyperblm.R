hyperblm <- function(formula, data, subset, weights, na.action,
                     xx = FALSE, y = FALSE, contrasts = NULL,
                     offset, method = "Nelder-Mead",
                     startMethod = "Nelder-Mead", startStarts = "BN",
                     paramStart = NULL,
                     maxiter = 100, tolerance = 0.0001,
                     controlBFGS = list(maxit = 1000),
                     controlNM = list(maxit = 10000),
                     maxitNLM = 10000,
                     controlCO = list(), silent = TRUE, ...){
    ret.x <- xx
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "offset"), names(mf), 0)
    mf <- mf[c(1,m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- model.weights(mf)
    offset <- model.offset(mf)
    if(!is.null(offset) && length(offset) != NROW(y))
        stop("Number of offsets is ", length(offset),
             ", should equal ", nrow(y), " (number of observations)")
    if(is.empty.model(mt)){
        x <- NULL
        regressionResult <- list(coefficients = numeric(0),
                                 distributionParams = numeric(0),
                                 residuals = y, fitted.values = 0 * y,
                                 weights = w, rank = 0,
                                 df.residual = length(y))
        if(!is.null(offset))
            regressionResult$fitted.values <- offset
    } else {
        x <- model.matrix(mt, mf, contrasts)
        regressionResult <- hyperblmFit(x, y, offset = offset,
                                        method = method,
                                        startMethod = startMethod,
                                        startStarts = startStarts,
                                        paramStart = paramStart,
                                        maxiter = maxiter,
                                        tolerance = tolerance,
                                        controlBFGS = controlBFGS,
                                        controlNM = controlNM,
                                        maxitNLM = maxitNLM,
                                        controlCO = controlCO,
                                        silent = silent, ...)
    }
    class(regressionResult) <- "hyperblm"
    regressionResult$na.action <- attr(mf, "na.action")
    regressionResult$offset <- offset
    regressionResult$contrasts <- attr(x, "contrasts")
    regressionResult$xlevels <- .getXlevels(mt, mf)
    regressionResult$call <- cl
    regressionResult$terms <- mt
    if(ret.x)
        regressionResult$x <- x
    if(ret.y)
        regressionResult$y <- y
    regressionResult
}

print.hyperblm <- function(x,
                           digits = max(3, getOption("digits")-3), ...){
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    coefficients <- x$coefficients
    distributionParams <- x$distributionParams
    cat("\nData:     ", x$obsName, "\n")
    if(length(coefficients) == 0){
        cat("\nNo coefficient\n")
    } else {
        cat("Regression coefficient estimates:\n")
        print.default(format(coefficients, digits = digits),
                      print.gap = 2, quote = FALSE)
    }
    if(length(distributionParams) == 0){
        cat("\nNo parameter for hyperbolic error distribution\n")
    } else {
        names(distributionParams) <- c("mu", "delta", "alpha", "beta")
        cat("Distribution Parameter estimates:\n")
        print.default(format(distributionParams, digits = digits),
                      print.gap = 2, quote = FALSE)
    }
    cat("Method:            ", x$method, "\n")
    cat("Likelihood:        ", x$mle, "\n")
    cat("Convergence code:  ", x$conv, "\n")
    cat("Iteration:         ", x$iter, "\n")
    invisible(x)
}

coef.hyperblm <- function(object, ...) {
    list(coefficients = object$coefficient,
         distributionParams = object$distributionParams)
}
