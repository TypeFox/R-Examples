## =============================================================================
## ----------------------------- declare generic -------------------------------
## =============================================================================
#' Extract Standard Errors of Model Coefficients
#'
#' These functions extract standard errors of model coefficients
#' from objects returned by count-modeling functions.
#'
#' When bootsrap standard error are required, the function cheks for the
#' bootsrap sample in \code{object}. If it is not found, the bootsrap
#' sample is created and a warning is sent.
#' @param object object returned by one of the count-modeling functions
#' @param parm parameters name or index
#' @param type type of standard error. User can choose between asymtotic
#' normal standard errors (\code{asymptotic}) or bootsrap (\code{boot}).
#' @param ... TODO
#' @return a named numeric vector
#' @export
se.coef <- function(object, parm, type, ...) {
    UseMethod("se.coef", object)
}

## ==============================================================================
## ------------------------------ methods definition ----------------------------
## ==============================================================================

#' Methods for renewal objects
#'
#' Methods for renewal objects.
#'
#' @param object,...,type,parm,level,bootType,x,digits TODO
#' @examples
#' fn <- system.file("extdata", "McShane_Wei_results_boot.RDS", package = "Countr")
#' object <- readRDS(fn)
#' class(object) # "renewal"
#'
#' coef(object)
#' vcov(object)
#'
#' ## Pearson residuals: rescaled by sd
#' head(residuals(object, "pearson"))
#' ## response residuals: not rescaled
#' head(residuals(object, "response"))
#'
#' head(fitted(object))
#'
#' ## loglik, nobs, AIC, BIC
#' c(loglik = as.numeric(logLik(object)), nobs = nobs(object),
#'   AIC = AIC(object), BIC = BIC(object))
#'
#' asym <- se.coef(object, , "asymptotic")
#' boot <- se.coef(object, , "boot")
#' cbind(asym, boot)
#' @name renewal_methods
NULL

#' @rdname renewal_methods
#' @method coef renewal
#' @export
coef.renewal <- function(object, ...) {
    object$coefficients
}

#' @rdname renewal_methods
#' @method vcov renewal
#' @export
vcov.renewal <- function(object, ...) {
    object$vcov
}

#' @rdname renewal_methods
#' @method residuals renewal
#' @export
residuals.renewal <- function (object, type = c("pearson", "response"), ...) {
    type <- match.arg(type)
    res <- object$residuals
    switch(type,
           "response" = {
               return(res)
           },
           "pearson" = {
               return(res / sqrt(object$wi))
           }
           )
}

#' @rdname renewal_methods
#' @method fitted renewal
#' @export
fitted.renewal <- function (object, ...) {
    object$fitted.values
}

#' Extract Standard Errors of Model Coefficients
#'
#' These functions extract standard errors of model coefficients
#' from objects returned by \code{renewal} functions.
#'
#' When bootsrap standard error are required, the function cheks for the
#' bootsrap sample in \code{object}. If it is not found, the bootsrap
#' sample is created and a warning is sent.
#' @param object object returned by \code{renewal}.
#' @param parm parameters name or index
#' @param type type of standard error. User can choose between asymtotic
#' normal standard errors (\code{asymptotic}) or bootsrap (\code{boot}).
#' @param ... TODO
#' @return an named numeric vector
#' @examples
#' ## see examples for renewal_methods
#' @export
se.coef.renewal <- function(object, parm, type = c("asymptotic", "boot"), ...) {
    type <- match.arg(type)
    switch(type,
           "asymptotic" = {
               se <- sqrt(diag(object$vcov))
               names(se) <- names(object$coefficients)
           },
           "boot" = {
               if (is.null(object$boot)) {
                   warning("boot sample not found. it will be created ...")
                   R <- list(...)$R
                   R <- ifelse(is.null(R), 100, R)
                   object <- addBootSampleObject.renewal(object, R = R)
               }
               ## extract se from car::summary.boot
               se <- summary(object$boot)[, "bootSE"]
           }
           )
    names(se) <- names(coef(object))

    ind <- seq(along = se)
    if (!missing(parm))
        ind <- which(parm %in% names(se))

    return(se[parm])
}


#' @method confint renewal
#' @examples
#' ## CI for coefficients
#' asym <- confint(object, type = "asymptotic")
#' ## Commenting out for now, see the nite in the code of confint.renewal():
#' ## boot <- confint(object, type = "boot", bootType = "norm")
#' ## list(asym = asym, boot = boot)
#' @rdname renewal_methods
#' @export
confint.renewal <- function(object, parm, level = 0.95,
                            type = c("asymptotic", "boot"),
                            bootType = c("norm", "bca", "basic", "perc", "all"),
                            ...) {
    type <- match.arg(type)
    bootType <- match.arg(bootType)
    switch(type,
           "asymptotic" = {
               return(confint.default(object, parm, level, ...))
           },
           "boot" = {
               if (is.null(object$boot)) {
                   warning("boot sample not found. it will be created ...")
                   R <- list(...)$R
                   R <- ifelse(is.null(R), 100, R)
                   object <- addBootSampleObject.renewal(object, R = R)
               }
               ## extract se from car::summary.boot
               ## 2016-03-23 TODO: something seems unfinished here
               ci <- confint(object$boot, parm, level, bootType, ...)
               rownames(ci) <- names(coef(object))
               return(ci)
           }
           )
}

#' @method summary renewal
#' @examples
#' summary(object)
#' @rdname renewal_methods
#' @export
summary.renewal <- function(object, ...) {
    object$residuals <- residuals(object, type = "pearson")
    se <- sqrt(diag(object$vcov))
    coef <- object$coefficients

    zstat <- coef/se
    pval <- 2 * pnorm(-abs(zstat))
    coef <- cbind(coef, se, zstat, pval)
    colnames(coef) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    object$coefficients <- coef
    object$fitted.values <- object$model <- object$y <- object$x <- NULL
    object$start <- NULL
    class(object) <- "summary.renewal"
    object
}

#' @method print renewal
#' @examples
#' print(object)
#' @rdname renewal_methods
#' @export
print.renewal <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") *
        0.85)), "", sep = "\n")
    if (!x$converged) {
        cat("model did not converge\n")
    } else {
        ##--------------- prepare links char
        textLink <- .summarizeLinkInformation(x$link)
        cat(paste0("\nCount model coefficients (inter-arrival ", x$dist,
                   " with ", textLink, "):\n"))
        print.default(format(x$coefficients, digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\nLog-likelihood:", sprintf(paste0("%.", digits, "f"), x$loglik),
            "on", x$n - x$df.residual, "Df\n")
    }
    invisible(x)
}

#' @rdname renewal_methods
#' @method print summary.renewal
#' @export
print.summary.renewal <- function(x, digits = max(3, getOption("digits") - 3),
                                  ...) {
    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") *
        0.85)), "", sep = "\n")
    if (!x$converged) {
        cat("model did not converge\n")
    } else {
        ## residuals print
        cat("Pearson residuals:\n")
        print(structure(quantile(x$residuals), names = c("Min",
            "1Q", "Median", "3Q", "Max")), digits = digits, ...)

        ## dist & link functions
        ##--------------- prepare links char
        textLink <- .summarizeLinkInformation(x$link)
        cat(paste0("\nCount model coefficients (inter-arrival ", x$dist,
                   " with ", textLink, "):\n"))

        ##---------------- print coefficients
        printCoefmat(x$coefficients, digits = digits, signif.legend = FALSE)
        if (getOption("show.signif.stars") &
            any(rbind(x$coefficients)[, 4] < 0.1)
            )
            cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
                "\n")

        ##------------------ optimization routine
        cat(paste("\nNumber of iterations in", x$method, "optimization:",
                  x$iterations, "\n"))

        ##-------------------- exec time
        cat(paste("\nExecution time", x$execTime, "\n"))

        ##------------------- log-likelihood
        cat("Log-likelihood:", sprintf(paste0("%.", digits, "f"), x$loglik),
            "on", x$n - x$df.residual, "Df\n")
    }
    invisible(x)
}

#' @rdname renewal_methods
#' @method model.matrix renewal
#' @export
model.matrix.renewal <- function(object, ...) {
    if (!is.null(object$x))
        return(object$x)
    else if (!is.null(object$model))
        return(.getModelMatrix(Formula(object$formula), object$dist,
                               object$model, object$anc)
               )
    else
        stop("model should be saved to create model.matrix !")
}

## \code{logLik} method for class \code{renewal}
## @param object an object from class \code{renewal}
## @param ... not used
#' @method logLik renewal
#' @rdname renewal_methods
#' @export
logLik.renewal <- function (object, ...) {
    structure(object$loglik, df = object$n - object$df.residual,
              class = "logLik")
}

## \code{nobs} method for class \code{renewal}
## @param object an object from class \code{renewal}
## @param ... not used
#' @method nobs renewal
#' @rdname renewal_methods
#' @examples
#' ## see renewal_methods
#' @export
nobs.renewal <- function (object, ...) {
    object$n
}

#' extractAIC method for class renewal
#'
#' extractAIC method for class renewal
#'
#' @param fit an object from class \code{renewal}
#' @param scale TODO
#' @param k TODO
#' @param ... not used
#' @method extractAIC renewal
#' @examples
#' ## see renewal_methods
#' @export
extractAIC.renewal <- function (fit, scale, k = 2, ...) {
     c(attr(logLik(fit), "df"), AIC(fit, k = k))
}

#' Predict method for renewal objects
#'
#' Compute predictions from renewal objects.
#'
#' @param type type of prediction.  If equal to \code{"response"}, give the mean
#'     probability associated with the individual covariates. If \code{"prob"},
#'     give the probability of the observed count.
#' @param time TODO
#' @inheritParams stats::predict.lm
#' @method predict renewal
#' @examples
#' fn <- system.file("extdata", "McShane_Wei_results_boot.RDS", package = "Countr")
#' object <- readRDS(fn)
#' data <- object$data
#' ## old data
#' predOld.response <- predict(object, type = "response", se.fit = TRUE)
#' predOld.prob <- predict(object, type = "prob", se.fit = TRUE)
#'
#' ## newData (extracted from old Data)
#' newData <- head(data)
#' predNew.response <- predict(object, newdata = newData,
#'                             type = "response", se.fit = TRUE)
#' predNew.prob <- predict(object, newdata = newData,
#'                         type = "prob", se.fit = TRUE)
#'
#' cbind(head(predOld.response$values),
#'            head(predOld.response$se$scale),
#'            head(predOld.response$se$shape),
#'            predNew.response$values,
#'            predNew.response$se$scale,
#'            predNew.response$se$shape)
#'
#' cbind(head(predOld.prob$values),
#'       head(predOld.prob$se$scale),
#'       head(predOld.prob$se$shape),
#'       predNew.prob$values,
#'       predNew.prob$se$scale,
#'       predNew.prob$se$shape)
#' @export
predict.renewal <- function(object, newdata = NULL, type = c("response", "prob"),
                            se.fit = FALSE, terms = NULL, na.action = na.pass,
                            time = 1.0, ...) {

    type <- match.arg(type)
    dist <- object$dist

    ## custom params
    customPars <- object$customPars

    ## link list
    linkList <- object$link

    ## check convolution parameters
    convPars <- renewal.convPars(list(...)$convPars, object$dist)

    ## prepare the modelMatrixList object as well as Y (response)
    if (is.null(newdata)) { ## no data provided
        ## build the modelMatrixList object
        if (!is.null(object$x))
            modelMatrixList <- object$x ## mf
        else if (!is.null(object$model)) {
            modelMatrixList <- .getModelMatrix(Formula(object$formula),
                                               object$dist,
                                               object$model, object$anc,
                                               customPars)
        } else
            stop(paste("predicted probabilities cannot be",
                       "computed with missing newdata")
                 )

        if (type == "prob") {
            ## check that the response is found in the object
            Y <- object$y
            if (is.null(Y))
                stop(paste("response should be saved in the fit object for",
                           "probability predictions !"))
        } else {
            out <- object$fitted.values
            se <- NA
            if (se.fit) {
                ## C <- .modelData(modelMatrixList, object$dist, customPars)
                ## se <- sqrt(diag(C %*% vcov(object) %*% t(C)))
                ## se_out <- .transformSE(se, object$link)
                se_out <- .getPredictionStd(modelMatrixList, vcov(object),
                                            object$dist, object$link, customPars)
            }

            return(list(values = as.numeric(out), se = se_out))
        }
    } else { ## data provided
        Fform <- Formula(object$formula)
        mf <- model.frame(Fform, data = newdata, na.action = na.action)
        modelMatrixList <- .getModelMatrix(Fform, object$dist, mf,
                                           object$anc, customPars)
        Y <- model.response(mf)
    }

    key <- ifelse(type == "response", TRUE, FALSE)

    out <-  .objectiveFunction(coef(object), dist, modelMatrixList,
                               linkList, time, convPars, Y, NULL, Ev = key,
                               seriesPars = object$seriesPars,
                               weiMethod = object$weiMethod,
                               summa = FALSE, customPars)

    if (type == "response") ## extract the mean response
        out <- sapply(out, .extractElem, ind = "ExpectedValue")

    se <- NA
    if (se.fit) {
        ## C <- .modelData(modelMatrixList, object$dist, customPars)
        ## se <- sqrt(diag(C %*% vcov(object) %*% t(C)))
        ## se_out <- .transformSE(se, object$link)
        ## se_out <- .getPredictionStd(modelMatrixList, vcov(object),
        ##                             object$dist, object$link, customPars)
        se_out <- .getPredictionStd(modelMatrixList, vcov(object),
                                    object$dist, object$link, customPars)
    }

    return(list(values = as.numeric(out), se = se_out))
}

#' Create a bootsrap sample for coefficient estimates
#'
#' Create a boostrap sample from coefficient estimates.
#'
#' The information in \code{object} is used to prepare the arguments and then
#' \code{boot} is called to generate the bootstrap sample.
#' The bootstrap sample is stored in \code{object} as component \code{"boot"}.
#'
#' @inheritParams se.coef
#' @param R numeric, number of bootsrap samples to generate.
#' @param ... not used.
#' @return \code{object} with additional component \code{"boot"}
#' @examples
#' ## see renewal_methods
#' @export
#' @importFrom boot boot
addBootSampleObject.renewal <- function(object, R = 100, ...) {
    formula <- object$formula
    data <- object$data
    weights <- object$weights
    dist <- object$dist
    anc <- object$anc
    convPars <- object$convPars
    link <- object$link
    time <- object$time
    control <- object$control
    control$trace <- FALSE
    customPars <- object$customPars

    bootIter <- 0

    bootFun <- function(data, indices) {
        bootIter <<- bootIter + 1
        cat("####### bootIter is: ", bootIter, "\n")

        bdata <- data[indices, , drop = FALSE]
        res <- renewal(formula = formula, data = bdata, weights = weights,
                       dist = dist, anc = anc, convPars = convPars,
                       link = link, time = time, computeHessian = FALSE,
                       control = control, customPars = customPars)
        coef(res)
    }

    object$boot <- boot(data = data, statistic = bootFun, R = R)
    object
}

#' @rdname renewal_methods
#' @method df.residual renewal
#' @export
df.residual.renewal <- function(object, ...) {
    object$df.residual
}

## copied update.default and modified it.
update.renewal <- function (object, formula., anc, ..., evaluate = TRUE)
{
    if (is.null(call <- getCall(object)))
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.))
        call$formula <- update.formula(formula(object), formula.)
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }

    ## process argument `anc' - this is the chunk specific to 'renewal'
    if(!missing(anc)  &&  !is.null(anc)) {
        newAnc <- call$anc # TODO: check if need as.list or whatever
        if(is.null(newAnc)) {
            call$anc <- anc
        } else {
            existingAncNames <- names(newAnc)
            ## TODO: for speed use match()
            for(nam in names(anc)) {
                if(nam %in% existingAncNames) {
                    newAnc[[nam]] <- update.formula(newAnc[[nam]], anc[[nam]])
                } else {
                    newAnc[[nam]] <- anc[[nam]]
                }
            }
            call$anc <- newAnc
        }
    }

    if (evaluate)
        eval(call, parent.frame())
    else call
}
