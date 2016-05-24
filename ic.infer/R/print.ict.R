print.ict <- function (x, digits = max(3, getOption("digits") - 3), scientific=FALSE, ...) 
{
    cat("\nOrder-related hypothesis test:\n")
        hilf <- c(x$TP, x$T, x$p.value)
        aus <- c(format(hilf[1], scientific=FALSE),
            format(hilf[2], digits = digits),
            format(if (hilf[3] < 0.0001) "<0.0001" else round(hilf[3], 4), digits=4))
        names(aus) <- c("Test problem", "Test statistic", "p-value")
        print.default(aus, print.gap = 2, quote = FALSE, ...)
    cat("\n")
    invisible(x)
}
