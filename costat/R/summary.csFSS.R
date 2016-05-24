summary.csFSS <-
function (object, size = 0.05, ...) 
{
    Nsims <- nrow(object$startpar)
    Ncoefs <- ncol(object$startpar)/2
    cat("Name of X time series: ", object$tsxname, "\n")
    cat("Name of Y time series: ", object$tsyname, "\n")
    cat("Length of input series: ", length(object$tsx), "\n")
    cat("There are ", Nsims, " sets of solutions\n")
    cat("Each solution vector is based on ", Ncoefs, " coefficients\n")
    if (all(object$convergence == 0)) 
        cat("Every solution converged o.k.\n")
    else cat("Some solutions did not converge, check convergence component for more information. Zero indicates successful convergence, other values mean different things and you should consult the help page for `optim' to discover what they mean\n")
    cat("For size level: ", size, "\n")
    cat("\t", tmp <- sum(object$pvals < size), " solutions appear NOT to be stationary\n")
    cat("\t", Nsims - tmp, " solutions appear to be stationary\n")
    r <- range(object$pvals)
    cat("Range of p-values: (", r[1], ",", r[2], ")\n")
    cat("\nWavelet filter for combinations: ", object$filter.number, 
        " ", object$family, "\n")
    cat("Wavelet filter for spectrum: ", object$spec.filter.number, 
        " ", object$spec.family, "\n")
}
