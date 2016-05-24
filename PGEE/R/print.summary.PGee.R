print.summary.PGee <-
function(x, digits = NULL, quote = FALSE, prefix = "", ... )
{
    if(is.null(digits))
        digits <- options()$digits
    else options(digits = digits)
    cat("\n",x$title)
    cat("\n",x$version,"\n")
    cat("\nModel:\n")
    cat(" Link:                     ",x$model$link,"\n")
    cat(" Variance to Mean Relation:",x$model$varfun,"\n")
    if(!is.null(x$model$M))
        cat(" Correlation Structure:    ",x$model$corstr,", M =",x$model$M,"\n")
    else 	cat(" Correlation Structure:    ",x$model$corstr,"\n")
    cat("\nCall:\n")
    dput(x$call)
    cat("\nSummary of Residuals:\n")
    print(x$residual.summary, digits = digits)
    nas <- x$nas
###	if(any(nas))
    if(!is.null(nas) && any(nas))
        cat("\n\nCoefficients: (", sum(nas),
            " not defined because of singularities)\n", sep = "")
    else cat("\n\nCoefficients:\n")
    print(x$coefficients, digits = digits)
    cat("\nEstimated Scale Parameter: ", format(round(x$scale,digits)))
    cat("\nLambda value: ", format(round(x$lambda.value, digits)))
    cat("\nNumber of Iterations: ", x$iterations)
    cat("\n\nWorking Correlation\n")
    print(x$working.correlation,digits=digits)
    if(!is.null(x$naive.correlation)) {
        cat("\n\nNaive Correlation of Estimates:\n")
        print(x$naive.correlation)
    }
    if(!is.null(x$robust.correlation)) {
        cat("\n\nRobust Correlation of Estimates:\n")
        print(x$robust.correlation)
    }
    invisible(x)
}
