##
## Analogue of fitAdaptiveLasso, for SCAD
##
##' @import ncvreg
fitSCAD <- function(X,
                    y,
                    polyTerms,
                    family,
                    lambda,
                    nlambda,
                    lambda.min.ratio,
                    nfolds,
                    thresh,
                    maxit,
                    ...)
{
    ## Compute polynomial expansion
    X.expand <- expandMatrix(X, polyTerms, intercept = FALSE)

    ## As in the cross-validation case of fitAdaptiveLasso(), need to use
    ## do.call() here to get around the fact that ncvreg() only checks for a
    ## missing `lambda` argument, not a NULL one
    ans <- list(X = X.expand,
                y = y,
                family = family,
                penalty = "SCAD",
                lambda.min = lambda.min.ratio,
                nlambda = nlambda,
                eps = thresh,
                max.iter = maxit,
                nfolds = nfolds)
    if (!is.null(lambda))
        ans$lambda <- lambda

    if (is.null(lambda) || length(lambda) > 1) {
        ## Cross-validate
        ans <- do.call(cv.ncvreg, ans)

        ## Extract information about cross-validation
        lambdaCV <- list(lambda = ans$lambda,
                         cvError = ans$cve,
                         lambdaMin = ans$lambda.min,
                         errorMin = min(ans$cve))

        ans <- list(coef = coef(ans),
                    lambda = lambdaCV$lambdaMin,
                    lambda.cv = lambdaCV)
                    
    } else {
        ## Fit directly using the specified penalization factor
        ans <- do.call(ncvreg, ans)
        ans <- list(coef = coef(ans),
                    lambda = lambda,
                    lambda.cv = NULL)
    }

    ans
}
