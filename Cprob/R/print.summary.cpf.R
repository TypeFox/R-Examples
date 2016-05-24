print.summary.cpf <- function(x, ...) {
    if (!inherits(x, "summary.cpf"))
        stop("'x' must be of class 'summary.cpf'")
    cat("Call: "); dput(x$call); cat("\n")
    if (!is.null(x$z)) {
        cat("Two sample test: \n")
        cat("\tz: ", x$z, "\n")
        cat("\tp: ", x$p, "\n\n")
    }
    for (i in seq_len(length(x$est))) {
        cat("\t", names(x$est[i]), "\n")
        ind <- floor(quantile(1:nrow(x$est[[i]]), prob = c(0, 0.25, 0.5, 0.75, 1)))
        print(x$est[[i]][ind, ], rownames = FALSE)
        cat("\n")
    }
    invisible(x)
}
