### Summary function for etmCIF objects

summary.etmCIF <- function(object, ci.fun = "cloglog", level = 0.95, ...) {

    if (!inherits(object, "etmCIF")) {
        stop("'object' must be of class 'etmCIF'")
    }

    l.X <- ncol(object$X)
    l.trans <- nrow(object[[1]]$trans)

    temp <- lapply(object[seq_len(l.X)], function(ll) {
        aa <- summary(ll, ci.fun = ci.fun, level = level, ...)[seq_len(l.trans)]
        names(aa) <- paste("CIF ", sapply(strsplit(sub("\\s", "|", names(aa)[1:l.trans]), "\\|"),
                                              "[", 2), sep = "")
        aa
    })

    class(temp) <- "summary.etmCIF"
    temp
}


### ... and the print function
print.summary.etmCIF <- function(x, ...) {

    if (!inherits(x, "summary.etmCIF")) {
        stop("'x' must be of class 'summary.etmCIF'")
    }

    for (i in seq_along(x)) {
        cat("\n\t", names(x)[i], "\n\n")
        time <- x[[i]][[1]]$time
        qtime <- quantile(time, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
        ind <- findInterval(qtime, time)
        for (j in seq_along(x[[i]])) {
            cat(names(x[[i]][j]), "\n")
            print(x[[i]][[j]][ind, ], row.names = FALSE)
            cat("\n")
        }
    }
    
    invisible()
}

        
