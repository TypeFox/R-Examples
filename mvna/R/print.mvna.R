print.mvna <- function(x, ...) {

    if (!inherits(x, "mvna")) {
        stop("'x' must be of class 'mvna'")
    }

    pos <- paste(x$trans[, 1], x$trans[, 2])
    temp <- x[pos]
    namen <- strsplit(pos, split = " ")

    for (i in seq_along(temp)) {
        time <- temp[[i]]$time
        qtime <- quantile(time, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
        ind <- findInterval(qtime, time, all.inside = TRUE)
        cat(paste("Transition", namen[[i]][1], "->", namen[[i]][2], "\n", sep = " "))
        print(round(temp[[i]][ind, ], 2), row.names = FALSE)
        cat("\n")
    }
    
    invisible()
    
}
