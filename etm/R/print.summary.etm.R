print.summary.etm <- function(x, ...) {
    if (!inherits(x, "summary.etm"))
        stop("'x' must be of class 'summary.etm'")
    if ("t" %in% names(x)) {
        cat(paste("No events between", x$s, "and", x$t, "\n\n", sep = " "))
        print(x$P[, , 1])
    }
    else {
        time <- x[[1]]$time
        qtime <- quantile(time, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
        ind <- findInterval(qtime, time)
            
        for (i in seq_along(x)) {
            cat(paste("Transition", names(x)[i], "\n", sep = " "))
            print(x[[i]][ind, ], row.names = FALSE)
            cat("\n")
        }
    }
    invisible()
}
