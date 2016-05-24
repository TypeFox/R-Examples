print.LORgee <-
function (x, ...) 
{
    cat("Call:\n")
    print(x$call)
    cat("\nNumber of clusters:", x$max.id, "\n")
    cat("Maximum cluster size:", max(x$clusz), "\n")
    cat("\nNumber of categories:", max(x$categories), "\n")
    cat("Number of observations:", x$nobs, "\n")
    cat("\nCoefficients:\n")
    print(x$coefficients)
    cat("\nNumber of iterations:", x$convergence$niter, "\n")
    cat("Algorithm converged:", x$convergence$conv, "\n")
    if (6 <= nrow(x$local.odds.ratios$theta)) {
        cat("\nLocal Odds Ratios Structure[1:6,1:6]\n")
        print(x$local.odds.ratios$theta[1:6, 1:6])
    }
    else {
        cat("\nLocal Odds Ratios Estimates[1:6,1:6]\n")
        print(x$local.odds.ratios$theta)
    }
    if (!is.null(x$pvalue)) 
        cat("\npvalue of Null model:", ifelse(x$pvalue<0.0001,"<0.0001", round(x$pvalue,4)), "\n")
}

