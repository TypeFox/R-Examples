#' @title print method for funreg object
#' @description Prints information from an object of class \code{funreg}.
#' @param x Object of class \code{funreg}.
#' @param digits Number of digits past the decimal place 
#' to show when printing quantities.
#' @param show.fits Whether to also print a table of fitted values for
#' the functional coefficients along a grid of times. 
#' @param ... Any other optional arguments to be passed to the default print
#' function.
#'@export
#'@method print funreg
print.funreg <- function(x, digits=4, show.fits=FALSE, ...) {
    stopifnot(class(x)=="funreg");
    summary1 <- summary(x,digits=digits);
    cat("funreg Functional Regression\n\n");
    cat("Call:\n");
    print(x$call.info);
    cat("\n");
    cat("Intercept estimate:      ");
    cat(x$intercept.estimate.uncentered);
    cat("\n");
    if (!is.null(x$other.covariates.estimate)) {
        cat("Subject-level coefficients:\n");
        print(summary1$subject.level.covariates.table,...);
    }
    if (show.fits) {
        cat("Functional coefficients:\n");
        print(summary1$functional.covariates.table,...);
    }
    if (!show.fits) {
        cat("\nTo view functional coefficients, use print(...,show.fits=TRUE)");
        cat(" or plot(...).\n\n");
    }
}