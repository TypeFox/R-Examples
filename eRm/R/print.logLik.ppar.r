`print.logLik.ppar` <-
function (x, digits = getOption("digits"),...)
{
    cat("Unconditional (joint) log Lik.: ", format(x$loglik, digits = digits),
        " (df=", format(x$df), ")\n", sep = "")
    invisible(x)
}
