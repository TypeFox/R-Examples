print.summary.mvna <- function(x, ...) {

    if (!inherits(x, "summary.mvna"))
        stop("'x' must be of class 'summary.mvna'")

    namen <- strsplit(names(x), split = " ")
    
    for (i in seq_along(x)) {
        time <- x[[i]]$time
        qtime <- quantile(time, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
        ind <- findInterval(qtime, time, all.inside = TRUE)
        cat(paste("Transition", namen[[i]][1], "->", namen[[i]][2], "\n", sep = " "))
        print(round(x[[i]][ind, ], 2), row.names = FALSE)
        cat("\n")
    }
    
    invisible()
    
}
