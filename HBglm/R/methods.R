################################
## Methods for hbglm objects   |
################################
##    * coef                   |
##    * print                  |
##    * summary                |
##    * print.summary          |
##    * predict                | # source code in predict.R
##    * plot                   | # source code in plot.R
################################

# Returns mean of coefficients' mcmc samples

coef.hbglm <- function(object, nburnin = 0, ...)
{
    nburnin <- if (nburnin <= 0) 0
    if (nburnin >= object$nsamples) stop("nburnin >= num samples.")
    stats <- sample.stats(object$samples, nburn = nburnin)
    
    coeffs <- list(beta = matrix(stats$beta[ , 1], nrow = object$model$J))
    rownames(coeffs$beta) <- object$model$grp.labels
    colnames(coeffs$beta) <- object$model$rand.cov
    if (object$family$has.tau) {
        coeffs$tau <- as.vector(stats$tau[ , 1])
        names(coeffs$tau) <- object$model$grp.labels
    }
    if (!is.null(stats$alpha)) {
        coeffs$alpha <- as.vector(stats$alpha[ , 1])
        names(coeffs$alpha) <- object$model$fixed.cov
    }
    if (!is.null(stats$theta)) {
        coeffs$theta <- matrix(stats$theta[ , 1], nrow = object$model$L)
        coeffs$Sigma <- matrix(stats$Sigma[ , 1], nrow = object$model$K)
        rownames(coeffs$theta) <- object$model$upper.cov
        colnames(coeffs$theta) <- object$model$rand.cov
        rownames(coeffs$Sigma) <- object$model$rand.cov
        colnames(coeffs$Sigma) <- object$model$rand.cov
    }
    return(coeffs)
}

# print method
print.hbglm <- function(x, digits = max(3, getOption("digits") -2),
                        width = getOption("width"), nburnin = 0, ...)
{
    coeff <- coef(x, nburnin = nburnin)
    cat("\n2-level Hierarchical Bayesian Regression model estimation ..\n")
    cat("\nMean coeffcient values are:\nbeta (random eff)\n")
    print(coeff$beta)
    if (!is.null(coeff$alpha)) {
        cat("\nalpha (fixed eff)\n"); print(coeff$alpha)
    }
    if (!is.null(coeff$tau)) {
        cat("\ntau\n"); print(coeff$tau)
    }
    if (!is.null(coeff$theta)) {
        cat("\ntheta (prior coeff)\n"); print(coeff$theta)
        cat("\nSigma (prior Covariance)\n"); print(coeff$Sigma)
    }
}

# Summary method
summary.hbglm <- function(object, nburnin = 0, ...)
{
    nburnin <- if (nburnin <= 0) 0
    if (nburnin >= object$nsamples) stop("nburnin >= num samples.")
    var.names <- make.var.names(object$model, object$family)
    tab <- sample.stats(object$samples, nburn = nburnin, var.names = var.names)
    object$Coeftable <- tab
    object$nburn <- nburnin
    class(object) <- c("summary.hbglm", "hbglm")
    return(object)
}

print.summary.hbglm <- function(x, digits = max(3, getOption("digits") -2),
                                  width = getOption("width"), ...)
{
    cat("\n2-level Hierarchical Bayesian Regression model estimation ..\n")
    cat("Model has:\n")
    cat(paste("\t", x$model$J, "groups\n"))
    cat(paste("\t", x$model$K, "random effects\n"))
    cat(paste("\t", x$model$M, "fixed effects\n"))
    cat(paste("\t", x$model$L, "upper level (prior) coefficients\n"))
    cat(paste0("Number of MCMC samples drawn = ", x$nsamples, "\n"))
    cat(paste0("Number of burn-in samples = ", x$nburn, "\n"))
    cat("\nCoefficients:\nbeta (random effects)\n")
    printCoefmat(x$Coeftable$beta, digits = digits)
    if (!is.null(x$Coeftable$alpha)) {
        cat("\nalpha (fixed eff)\n"); printCoefmat(x$Coeftable$alpha)
    }
    if (!is.null(x$Coeftable$tau)) {
        cat("\ntau\n"); printCoefmat(x$Coeftable$tau)
    }
    if (!is.null(x$Coeftable$theta)) {
        cat("\ntheta (prior coeff)\n"); printCoefmat(x$Coeftable$theta)
        cat("\nSigma (prior Covariance)\n"); printCoefmat(x$Coeftable$Sigma)
    }
}
