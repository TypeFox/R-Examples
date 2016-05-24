summary.cv.grpreg <- function(object, ...) {
  S <- pmax(object$null.dev - object$cve, 0)
  rsq <- S/object$null.dev
  snr <- S/object$cve
  nvars <- predict(object$fit, type="nvars")
  ngroups <- predict(object$fit, type="ngroups")
  model <- switch(object$fit$family,
                  gaussian = "linear",
                  binomial = "logistic",
                  poisson = "Poisson")
  d <- dim(object$fit$beta)
  if (length(d)==3) {
    p <- d[2] - 1
  } else {
    p <- d[1] - 1
  }
  val <- list(penalty=object$fit$penalty, model=model, n=object$fit$n, p=p, min=object$min, lambda=object$lambda, cve=object$cve, r.squared=rsq, snr=snr, nvars=nvars, ngroups=ngroups, d=d)
  if (object$fit$family=="gaussian") val$sigma <- sqrt(object$cve)
  if (object$fit$family=="binomial") val$pe <- object$pe
  structure(val, class="summary.cv.grpreg")
}
print.summary.cv.grpreg <- function(x, digits, ...) {
  digits <- if (missing(digits)) digits <- c(2, 4, 2, 2, 3) else rep(digits, length.out=5)
  if (length(x$d)==3) {
    cat(x$penalty, "-penalized multivariate ", x$model, " regression with m=", x$d[1], ", n=", x$n/x$d[1], ", p=", x$p, "\n", sep="")
  } else {
    cat(x$penalty, "-penalized ", x$model, " regression with n=", x$n, ", p=", x$p, "\n", sep="")
  }
  cat("At minimum cross-validation error (lambda=", formatC(x$lambda[x$min], digits[2], format="f"), "):\n", sep="")
  cat("-------------------------------------------------\n")
  cat("  Nonzero coefficients: ", x$nvars[x$min], "\n", sep="")
  cat("  Nonzero groups: ", x$ngroups[x$min], "\n", sep="")
  cat("  Cross-validation error of ", formatC(min(x$cve), digits[1], format="f"), "\n", sep="")
  cat("  Maximum R-squared: ", formatC(max(x$r.squared), digits[3], format="f"), "\n", sep="")
  cat("  Maximum signal-to-noise ratio: ", formatC(max(x$snr), digits[4], format="f"), "\n", sep="")
  if (x$model == "logistic") cat("  Prediction error at lambda.min: ", formatC(x$pe[x$min], digits[5], format="f"), "\n", sep="")
  if (x$model == "linear") cat("  Scale estimate (sigma) at lambda.min: ", formatC(sqrt(x$cve[x$min]), digits[5], format="f"), "\n", sep="")
}
