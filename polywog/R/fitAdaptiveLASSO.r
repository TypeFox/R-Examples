##
## Calculate adaptive LASSO results.  k-fold cross-validation is used to select
## the penalization factor if `lambda` is NULL or contains multiple values.
##
##' @import glmnet
fitAdaptiveLASSO <- function(X,
                             y,
                             weights,
                             polyTerms,
                             family,
                             penwt,
                             lambda,
                             nlambda,
                             lambda.min.ratio,
                             nfolds,
                             foldid,
                             thresh,
                             maxit,
                             .parallel,
                             ...)
{
    ## Compute polynomial expansion
    X.expand <- expandMatrix(X, polyTerms, intercept = FALSE)

    if (is.null(lambda) || length(lambda) > 1) {
        ## Cross-validate

        ## cv.glmnet() will fail if called with `foldid = NULL`, since it only
        ## checks for *missing* `foldid` to construct fold IDs itself, so we
        ## need to call the function indirectly
        ans <- list(x = X.expand,
                    y = y,
                    family = family,
                    weights = weights,
                    nlambda = nlambda,
                    lambda.min.ratio = lambda.min.ratio,
                    lambda = lambda,
                    standardize = FALSE,
                    thresh = thresh,
                    penalty.factor = penwt,
                    maxit = maxit,
                    nfolds = nfolds,
                    parallel = .parallel)
        if (!is.null(foldid))
            ans$foldid <- foldid
        ans <- do.call(cv.glmnet, ans)

        ## Extract information about cross-validation
        lambdaCV <- list(lambda = ans$lambda,
                         cvError = ans$cvm,
                         lambdaMin = ans$lambda.min,
                         errorMin = min(ans$cvm))

        ans <- list(coef = coef(ans, s = "lambda.min")[, 1],
                    lambda = lambdaCV$lambdaMin,
                    lambda.cv = lambdaCV)
    } else {
        ## Fit directly using the specified penalization factor
        ans <- glmnet(x = X.expand,
                      y = y,
                      family = family,
                      weights = weights,
                      lambda = lambda,
                      standardize = FALSE,
                      thresh = thresh,
                      penalty.factor = penwt,
                      maxit = maxit)
        ans <- list(coef = coef(ans)[, 1],
                    lambda = lambda,
                    lambda.cv = NULL)
    }

    ans
}
