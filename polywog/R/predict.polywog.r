##' Predict method for polywog objects
##'
##' Generates fitted values, including bootstrap confidence intervals, for in-
##' and out-of-sample data from a fitted polywog model.
##' @param object a fitted model of class \code{"polywog"}, typically the output
##' of \code{\link{polywog}}.
##' @param newdata an optional data frame containing observations for which
##' fitted values should be computed.  If not specified, fitted values are
##' generated for the data used to fit the model.
##' @param type specifies whether the fitted values should be generated on the
##' link scale (\eqn{X \beta}) or in terms of the expected value of the response
##' variable.  These only differ for binomial family models.
##' @param interval logical: whether to calculate bootstrap confidence intervals
##' for each fitted value.
##' @param level confidence level for the intervals.
##' @param bag logical: whether to use "bootstrap aggregation" to generate the
##' main fitted values (if \code{FALSE}, they are calculated from the main model
##' fit).
##' @param na.action a function specifying what to do with observations in
##' \code{newdata} containing \code{NA}s (default \code{\link{na.pass}}).  See
##' "Details".
##' @param ... other arguments, currently ignored.
##' @return If \code{interval = TRUE}, a matrix containing each fitted value and
##' its confidence interval.  Otherwise, a vector containing the fitted values.
##' @seealso For more user-friendly generation of fitted values, see
##' \code{\link{predVals}}.  To compute marginal effects, see
##' \code{\link{margEff.polywog}}.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @method predict polywog
##' @export
predict.polywog <- function(object, newdata,
                            type = c("link", "response"),
                            interval = FALSE, level = .95,
                            bag = FALSE,
                            na.action = na.pass, ...)
{
    type <- match.arg(type)
    transform <- (type == "response" && object$family == "binomial")

    ## Check for nothing that depends on bootstrap results if 'object' does not
    ## have a 'boot.matrix' element
    if (is.null(object$boot.matrix) && (interval || bag))
    {
        interval <- bag <- FALSE
        warning("Options 'interval' and 'bag' not available for models without a 'boot.matrix' element")
    }

    ## Setup largely adapted from predict.lm() code; the bits relating to
    ## 'newdata' being a model frame are adapted from mgcv::predict.gam()
    X.exists <- FALSE
    if (missing(newdata) || is.null(newdata)) {
        ## Use original model matrix or model frame if available
        X <- object$X
        if (is.null(X)) {
            if (is.null(object$model))
                stop("Fitted object must contain either 'model' or 'X' to use predict.polywog without specifying 'newdata'; re-run polywog with \"model = TRUE\"")
            mf <- object$model
        } else {
            X.exists <- TRUE
        }
        nd.is.mf <- FALSE
    } else if (is.data.frame(newdata) && !is.null(attr(newdata, "terms"))) {
        ## 'newdata' is a model frame -- this case must be treated separately,
        ## or else predVals() and margEff.polywog() won't work when the
        ## original model formula contains transformations of the original
        ## inputs

        nd.is.mf <- TRUE
    } else {
        ## Construct model frame from 'newdata'
        Terms <- delete.response(terms(object))
        mf <- model.frame(Terms, newdata, na.action = na.action, xlev =
                          object$xlevels)

        ## Check validity of 'newdata' (covariate types same as in fitted
        ## model)
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, mf)

        nd.is.mf <- FALSE
    }

    ## Compute the model matrix
    if (!X.exists) {
        X <- makeX(object$formula, if (nd.is.mf) newdata else mf)
    }

    ## Call the C++ backend to compute the predicted values and (if requested)
    ## confidence intervals
    pred <- computePredict(X = X,
                           poly_terms = object$polyTerms,
                           coef = list(main = coef(object),
                           boot = if (interval || bag) object$boot.matrix),
                           forPredVals = FALSE,
                           interval = interval,
                           bag = bag,
                           level = level,
                           transform = transform)
    if (interval) {
        pred <- do.call(cbind, pred)
    } else {
        pred <- pred$fit
    }

    ## If just computing in-sample fits, ensure conformity with the original
    ## 'na.action' (i.e. padding when na.action = na.exclude, as in
    ## 'predict.lm')
    if (missing(newdata) || is.null(newdata)) {
        pred <- napredict(object$na.action, pred)
    } else if (nd.is.mf) {
        pred <- napredict(attr(newdata, "na.action"), pred)
    } else {
        pred <- napredict(attr(mf, "na.action"), pred)
    }

    return(pred)
}
