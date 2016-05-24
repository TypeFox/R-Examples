## function: summary cv.grpregOverlap, need revise
# ------------------------------------------------------------------------------
summary.cv.grpregOverlap <- function(object, ...) {
    obj.new <- object
    class(obj.new) <- 'cv.grpreg'
    res <- summary(obj.new, ...)
    
    nvars.latent <- predict(object$fit, type="nvars", latent = T)
    res$nvars.latent <- nvars.latent
    d.latent <- dim(object$fit$beta.latent)
    if (length(d.latent)==3) {
        p.latent <- d.latent[2] - 1
    } else {
        p.latent <- d.latent[1] - 1
    }
    res$p.latent <- p.latent
    class(res) <- c('summary.cv.grpregOverlap', 'summary.cv.grpreg')
    res
}
# -------------------------------------------------------------------------------

## function: print summary.cv.grpregOverlap
# -------------------------------------------------------------------------------
print.summary.cv.grpregOverlap <- function(x, digits, ...) {
    if (missing(digits)) {
        digits <- c(2, 4, 2, 2, 3)
    } else {
        digits <- rep(digits, length.out=5)
    }
    cat("---------------------------------------------------------------\n")
    cat("Note: Overlapping-group selection via penalized regression: \n")
    cat("      'p.latent' is the number of latent variables!\n")
    cat("---------------------------------------------------------------\n")
    if (length(x$d)==3) {
        cat(x$penalty, "-penalized multivariate ", x$model, 
            " regression with m=", x$d[1], ", n=", x$n/x$d[1], ", p=", x$p, 
            " (p.latent=", x$p.latent, ")", "\n", sep="")
    } else {
        cat(x$penalty, "-penalized ", x$model, " regression with n=", x$n, 
            ", p=", x$p, " (p.latent=", x$p.latent, ")", "\n", sep="")
    }
    cat("At minimum cross-validation error (lambda=", formatC(x$lambda[x$min], digits[2], format="f"), "):\n", sep="")
    cat("---------------------------------------------------------------\n")
    cat("  Nonzero        coefficients: ", x$nvars[x$min], "\n", sep="")
    cat("  Nonzero latent coefficients: ", x$nvars.latent[x$min], "\n", sep="")
    cat("  Nonzero groups: ", x$ngroups[x$min], "\n", sep="")
    cat("  Cross-validation error of ", 
        formatC(min(x$cve), digits[1], format="f"), "\n", sep="")
    cat("  Maximum R-squared: ", 
        formatC(max(x$r.squared), digits[3], format="f"), "\n", sep="")
    cat("  Maximum signal-to-noise ratio: ", 
        formatC(max(x$snr), digits[4], format="f"), "\n", sep="")
    if (x$model == "logistic") {
        cat("  Prediction error at lambda.min: ", 
            formatC(x$pe[x$min], digits[5], format="f"), "\n", sep="")
    }
    if (x$model == "linear") {
        cat("  Scale estimate (sigma) at lambda.min: ", 
            formatC(sqrt(x$cve[x$min]), digits[5], format="f"), "\n", sep="")
    }
}
# -------------------------------------------------------------------------------
