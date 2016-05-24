print.cpf <- function(x, ...) {
    if (!inherits(x, "cpf"))
        stop("'x' must be of class 'cpf'")
    cat("Call: "); dput(x$call); cat("\n")
    cat("Number of observations: ",
        sum(x$n.risk[c(1, x$size.strata[-length(x$size.strata)] + 1)]),
        "\n\n")
    if (!is.null(x$X)) {
        cat("Covariate: ", colnames(x$X), "\n")
        cat("\tlevels: ", x$X[,], "\n\n")
    }
    if (!is.null(x$p)) {
        cat("Two sample test: \n")
        cat("\tz: ", x$z, "\n")
        cat("\tp: ", x$p, "\n\n")
    }
    Cause <- c(x$failcode, "other", "censored")
    n.event <- c(sum(x$n.event[, 1]), sum(x$n.event[, 2]), sum(x$n.lost))
    print(data.frame(Cause, n.event), row.names = FALSE)
    invisible()
}
