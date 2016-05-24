`print.logLik.eRm` <-
function (x, digits = getOption("digits"),...)
{
    cat("Conditional log Lik.: ", format(x$loglik, digits = digits),
        " (df=", format(x$df), ")\n", sep = "")
    invisible(x)
}

