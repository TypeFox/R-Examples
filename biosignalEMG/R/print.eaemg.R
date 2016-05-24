print.eaemg <- function(x, ...) {
    object <- x
    cat("Envelope averaged EMG")
    if (object$empirical) {
        cat("\n", round(object$level * 100, 2), "% empirical prediction intervals")
    } else {
        cat("\n", round(object$level * 100, 2), "% empirical prediction intervals")
    }
    cat("\n", dim(object$intervals)[1], "group points\n")
    cat("\nPrediction intervals:\n")
    print(object$intervals)
} 
