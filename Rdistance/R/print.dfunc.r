print.dfunc <- function( x, ... ){
#
#   Print a distance function
#

    cat("Call: ", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x)), print.gap = 2,
            quote = FALSE)
    }
    else cat("No coefficients\n")

    cat("\n")

    if( x$convergence == 0 ){
        mess <- "Success"
    } else {
        mess <- paste( "FAILURE (Exit code=", x$convergence, ", ", x$fit$message, ")")
    }
    cat(paste("Convergence: ", mess,  "\n", sep=""))


    if( x$expansions==0 ){
        mess <- ""
    } else {
        mess <- paste( "with", x$expansions, "expansion(s) of", casefold( x$series, upper=TRUE ), "series")
    }
    cat(paste("Function:", casefold(x$like.form, upper=TRUE), mess, "\n") )

    cat(paste("Strip:", x$w.lo, "to", x$w.hi, "\n"))
    cat(paste("Effective strip width:", format(ESW(x)), "\n"))
    cat(paste("Scaling: g(", x$x.scl, ") = ", format(x$g.x.scl), "\n", sep=""))
    cat(paste("Log likelihood:", format(x$loglik), "\n"))
    cat(paste("AIC:", format(AIC(x)), "\n"))


    cat("\n")
    invisible(x)
}
