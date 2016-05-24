`print.bootstrap.wa` <- function(x, ...) {
    s <- strsplit(x$deshrink, " ")[[1]]
    deshrink <- paste(toupper(substring(s, 1,1)), substring(s, 2),
                      sep="", collapse=" ")
    CV <- x$CV.method
    ## printing
    cat("\n")
    writeLines(strwrap("Weighted Averaging Bootstrap Predictions",
                       prefix = "\t"))
    cat("\nCall:\n")
    cat(paste(deparse(x$call), "\n\n"))
    cat(paste("Deshrinking     :", deshrink, "\n"))
    cat(paste("Crossvalidation :", CV, "\n"))
    cat(paste("Tolerance DW    :",
              ifelse(x$tol.dw, "Yes", "No"), "\n"))
    #cat(paste("No. samples     :", x$n.samp, "\n"))
    #cat(paste("No. species     :", x$n.spp, "\n\n"))
    perf <- performance(x)
    perf.names <- names(perf)
    attributes(perf) <- NULL
    names(perf) <- perf.names
    cat("\nPerformance:\n")
    print.default(round(perf, 4), print.gap = 2, ...)
    cat("\nTraining set predictions:\n")
    print.default(round(x$model$pred, 4), print.gap = 1, ...)
    invisible(x)
}
