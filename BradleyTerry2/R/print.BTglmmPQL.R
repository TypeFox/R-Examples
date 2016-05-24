print.BTglmmPQL <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    if (identical(x$sigma, 0)){
        cat("PQL algorithm converged to fixed effects model\n")
        return(NextMethod())
    }
    cat("\nCall: ", deparse(x$call), "\n", sep = "", fill = TRUE)
    if (length(coef(x))) {
        cat("Fixed effects:\n\n")
        print.default(format(x$coefficients, digits = digits),
            print.gap = 2, quote = FALSE)
    }
    else cat("No fixed effects\n\n")
    cat("\nRandom Effects Std. Dev.:", x$sigma, "\n")
    if (nzchar(mess <- naprint(x$na.action)))
        cat("\n", mess, "\n", sep = "")
}
