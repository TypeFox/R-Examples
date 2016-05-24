print.bic.glm <-
function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("\n Posterior probabilities(%): \n")
    out <- x$probne0
#   names(out) <- if (x$factor.type) x$namesx else colnames(x$mle)[-1]
    print(out, ...)
    cat("\n Coefficient posterior expected values: \n")
    out <- x$postmean
#    names(out) <- c("(Intercept)", colnames(x$mle)[-1])
    fout <- format(out, digits = digits)
    fout[is.na(out)] <- ""
    print.default(fout, print.gap = 2, quote = FALSE, ...)
    invisible(x)
}

