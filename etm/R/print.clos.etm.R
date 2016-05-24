print.clos.etm <- function(x, ...) {
    if (!inherits(x, "clos.etm")) {
        stop("'x' must be of class 'clos.etm'")
    }
    cat("The expected change in length of stay is:\n")
    cat(paste(round(x$e.phi, 3)), "\n")
    if (!is.null(x$e.phi.weights.1)) {
        cat("\nAlternative weighting:\n\n")
        cat(paste("Expected change in LOS with weight.1:",
                  round(x$e.phi.weights.1, 3), "\n", sep = " "))
        cat(paste("Expected change in LOS with weight.other:",
                  round(x$e.phi.weights.other, 3), "\n", sep = " "))
    }
    invisible()
}
