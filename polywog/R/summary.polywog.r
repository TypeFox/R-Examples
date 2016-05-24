##' Summarize a fitted polywog model
##'
##' Generates a "regression table" to summarize the fitted model, including
##' coefficients along with their bootstrapped standard errors and confidence
##' intervals.  If the fitted model does not have a \code{boot.matrix} element,
##' the output will contain \code{NA}s for the standard errors, and confidence
##' intervals will not be displayed.
##' @param object a fitted model of class \code{"polywog"}, typically the output
##' of \code{\link{polywog}}.
##' @param level width of the bootstrap confidence interval to compute for the
##' model coefficients.
##' @param prop0 logical: whether to print the proportion of bootstrap
##' iterations in which each coefficient was estimated as exactly 0.  This may
##' be informative but should \emph{not} be interpreted as a p-value.
##' @param ... other arguments, currently ignored.
##' @return An object of class \code{"summary.polywog"} whose elements are the
##' "regression table" (\code{coefficients}) and additional information from the
##' original fitted model.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @method summary polywog
##' @export
##' @importMethodsFrom Matrix rowMeans
summary.polywog <- function(object, level = .95, prop0 = FALSE, ...)
{
    ans <- list()

    ## Create matrix of summary results to be printed
    cf <- coef(object)
    se <- sqrt(diag(vcov(object)))
    ans$coefficients <- cbind("Estimate" = cf, "Std. Error" = se)
    if (!is.null(object$boot.matrix)) {
        q <- 0.5 - (level/2)
        interval <- apply(object$boot.matrix, 1, quantile, probs = c(q, 1-q))
        p0 <- rowMeans(object$boot.matrix == 0)
        ans$coefficients <- cbind(ans$coefficients, t(interval),
                                  "Prop. 0" = if (prop0) p0)
    }

    ## Add relevant information for printing model summary
    ans$call <- object$call
    ans$degree <- object$degree
    ans$family <- object$family
    ans$method <- object$method
    ans$penwt.method <- object$penwt.method
    ans$nobs <- object$nobs
    ans$nboot <-
        if (!is.null(object$boot.matrix)) ncol(object$boot.matrix) else 0
    ans$lambda <- object$lambda

    class(ans) <- "summary.polywog"
    return(ans)
}

##' @S3method print summary.polywog
print.summary.polywog <- function(x, digits = max(3, getOption("digits") - 3),
                                  ...)
{
    ## Print the function call used to fit the model (using same code as in
    ## 'print.lm')
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")

    ## Print a typical 'summary'-style matrix
    cat("Coefficients:\n")
    printCoefmat(coef(x), digits = digits, ...)

    ## Print important info about the model
    cat("\nRegularization method:",
        switch(x$method, alasso = "Adaptive LASSO", scad = "SCAD"))
    if (x$method == "alasso") {
        pname <- switch(x$penwt.method,
                        lm = "inverse linear model coefficients",
                        glm = "inverse logistic regression coefficients")
        cat("\nAdaptive weights:", pname)
    }
    cat("\nNumber of observations:", x$nobs)
    cat("\nPolynomial expansion degree:", x$degree)
    cat("\nModel family:", x$family)
    cat("\nBootstrap iterations:", x$nboot)
    cat("\nPenalization parameter (lambda):", format(x$lambda, digits = digits))
    cat("\n\n")

    invisible(x)
}
