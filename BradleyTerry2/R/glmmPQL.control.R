glmmPQL.control <- function (maxiter = 50, IWLSiter = 10, tol = 1e-6,
                             trace = FALSE) {
    call <- as.list(match.call())

    if (length(call) > 1) {
        argPos <- match(c("maxiter", "IWLSiter", "tol"), names(call))
        for (n in argPos[!is.na(argPos)]) {
            if (!is.numeric(call[[n]]) || call[[n]] <= 0)
                stop("value of '", names(call)[n], "' must be > 0")
        }
    }

    list(maxiter = maxiter, IWLSiter = IWLSiter, tol = tol,
         trace = trace)
}
