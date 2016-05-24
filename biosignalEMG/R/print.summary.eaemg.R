print.summary.eaemg <- function(x, ...) {
    object <- x
    cat("Envelope averaged EMG")
    if (object$empirical) {
        cat("\n", round(object$level * 100, 2), "% empirical prediction intervals")
    } else {
        cat("\n", round(object$level * 100, 2), "% empirical prediction intervals")
    }
    cat("\n", object$gp, "group points\n")
} 
