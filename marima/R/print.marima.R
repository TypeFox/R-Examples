##' @title print.marima
##' 
##' @description Print some (most relevant) content of a marima object.
##'
##' @param x = a marima object with results of marima analysis.
##' @param estimates = TRUE/FALSE: printout of parameter estimates.
##' @param pvalues   = TRUE/FALSE: printout of (approximate) p-values for
##' parameter estimates.
##' @param pattern   = TRUE/FALSE: printout of model definition pattern(s).
##' @param fvalues   = TRUE/FALSE: printout of parameter (approximate)
##' F-values.
##' @param ... Not used.
##'
##' @method print marima
##' @export

print.marima <- function(x, estimates = TRUE, pvalues = FALSE,
                         pattern = TRUE, fvalues = TRUE, ... ) {

    cat("Model dimension = kvar = ", x$kvar, " N = ", x$N, " \n")

    # if(pattern==TRUE){ cat('ar.pattern = \n')
    # print(short.form(x$call.ar.pattern,leading=F))
    # cat('ar.pattern = \n')
    # print(short.form(x$call.ma.pattern,leading=F))}

    cat("Averages for all variables:\n", x$averages, "\n")
    cat("Covariance(all data):\n")
    print(round(x$data.cov, 4))
    cat("Covariance(residuals):\n")
    print(round(x$resid.cov, 4))
    cat("Random variables are", x$randoms, "\n")

    if (pattern == TRUE) {
        cat("\n")
        cat("AR definition:\n")

        AR <- short.form(x$call.ar.pattern, leading = FALSE)
        print(AR)

        cat("MA definition:\n")
        MA <- short.form(x$call.ma.pattern, leading = FALSE)
        print(MA)

        if (pattern != TRUE) {
            cat("No model identification output specified \n")
        }
    }

    if (estimates == TRUE) {
        cat("AR estimates:\n")
        print(round(short.form(x$ar.estimates, leading = FALSE), 4))

        if (fvalues == TRUE) {
            cat("AR f-values (squared t-values):\n")
            print(round(short.form(x$ar.fvalues, leading = FALSE), 4))
        }
        if (pvalues == TRUE) {
            cat("AR p-values (%):\n")
            print(round(100 * short.form(x$ar.pvalues, leading = FALSE), 2))
        }
    }

    if (estimates != TRUE) {
        cat("No estimated model output specified \n")
    }

    if (estimates == TRUE) {
        cat("MA estimates:\n")
        print(round(short.form(x$ma.estimates, leading = FALSE), 4))
        if (fvalues == TRUE) {
            cat("MA f-values (squared t-values):\n")
            print(round(short.form(x$ma.fvalues, leading = FALSE), 2))
        }
        if (pvalues == TRUE) {
            cat("MA p-values (%):\n")
            print(round(100 * short.form(x$ma.pvalues, leading = FALSE), 2))
        }
    }
    if (estimates != TRUE) {
        cat("No estimated model output specified \n")
    }

}
