print.summary.mpp <- function(x, ...){
    cat("\nSummary of Marked Point Process Model\n\n")
    cat("Parameter Estimates:\n")
    cat(x$params, "\n\n")
    cat("Time Interval (TT):\n")
    cat(x$TT, "\n\n")
    cat("Ground Intensity Function Parameter Mapping (gmap):\n")
    print(x$gmap)
    cat("\n")
    cat("Mark Distribution Parameter Mapping (mmap):\n")
    print(x$mmap)
    cat("\n")
    cat("Log-Likelihood =", x$LL, "\n\n")
}

