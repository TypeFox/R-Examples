print.PGee <-
function(x, digits = NULL, quote = FALSE, prefix = "", ...)
{
    if(is.null(digits)) digits <- options()$digits else options(digits =
                                                                digits)
    cat("\n", x$title)
    cat("\n", x$version, "\n")
    cat("\nModel:\n")
    cat(" Link:                     ", x$model$link, "\n")
    cat(" Variance to Mean Relation:",x$model$varfun,"\n")
    if(!is.null(x$model$M))
        cat(" Correlation Structure:    ", x$model$corstr, ", M =", x$
            model$M, "\n")
    else cat(" Correlation Structure:    ", x$model$corstr, "\n")
    cat("\nCall:\n")
    dput(x$call)                        #       cat("\nTerms:\n")
###        ys <- matrix(rep(as.matrix(x$id, ncol = 1), 5), ncol = 5)
    ys <- matrix(rep(matrix(x$id, ncol = 1), 5), ncol = 5)
    ys[, 2] <- x$y
    ys[, 3] <- x$linear.predictors
    ys[, 4] <- x$fitted.values
    ys[, 5] <- x$residuals
    dimnames(ys) <- list(1:length(x$y), c("ID", "Y", "LP", "fitted",
                                          "Residual")) #       cat("\nFitted Values:\n")
    cat("\nNumber of observations : ", x$nobs, "\n")
    cat("\nMaximum cluster size   : ", x$max.id, "\n")
    nas <- x$nas
    if(any(nas))
        cat("\n\nCoefficients: (", sum(nas),
            " not defined because of singularities)\n", sep = "")
    else cat("\n\nCoefficients:\n")
    print(x$coefficients, digits = digits)
    cat("\nEstimated Scale Parameter: ", format(round(x$scale, digits)))
    cat("\nLambda value: ", format(round(x$lambda.value, digits)))
    cat("\nNumber of Iterations: ", x$iterations)
    cat("\n\nWorking Correlation[1:4,1:4]\n")
    print(x$working.correlation[1:4, 1:4], digits = digits)
    cat("\n\nReturned Error Value:\n")
    print(x$error)
    invisible(x)
}
