##' Experimental data from the job training program first studied by LaLonde (1986)
##'
##' A dataset of units in an experimental evaluation of a job training
##' program. Subset to those units with two years of pre-treatment income data.
##'
##'
##' \itemize{
##'   \item \code{age} - age in years.
##'   \item \code{education} - number of years of schooling.
##'   \item \code{black} - 1 if black, 0 otherwise.
##'   \item \code{hispanic} - 1 if Hispanic, 0 otherwise.
##'   \item \code{married} - 1 if married, 0 otherwise.
##'   \item \code{nodegree} - 1 if no high school degree, 0 otherwise.
##'   \item \code{re74} - earnings ($) in 1974.
##'   \item \code{re75} - earnings ($) in 1975.
##'   \item \code{re78} - earnings ($) in 1978.
##'   \item \code{u74} - 1 if unemployed in 1974, 0 otherwise.
##'   \item \code{u75} - 1 if unemployed in 1975, 0 otherwise.
##'   \item \code{treat} - 1 if treated, 0 otherwise.
##' }
##'
##' @docType data
##' @keywords datasets
##' @name lalonde.exp
##' @usage data(lalonde.exp)
##' @format A data frame with 445 rows and 12 variables
##' @references LaLonde, Robert J. (1986). Evaluating the Econometric
##' Evaluations of Training Programs with Experimental Data. The
##' American Economic Review, 76(4), 604--620.
NULL

##' Non-experimental data from Lalonde (1986)
##'
##' A dataset of experimental treated units and non-experimental
##' control units from the Panel Study of Income Dynamics (PSID).
##'
##' \itemize{
##'   \item \code{age} - age in years.
##'   \item \code{education} - number of years of schooling.
##'   \item \code{black} - 1 if black, 0 otherwise.
##'   \item \code{hispanic} - 1 if Hispanic, 0 otherwise.
##'   \item \code{married} - 1 if married, 0 otherwise.
##'   \item \code{nodegree} - 1 if no high school degree, 0 otherwise.
##'   \item \code{re74} - earnings ($) in 1974.
##'   \item \code{re75} - earnings ($) in 1975.
##'   \item \code{re78} - earnings ($) in 1978.
##'   \item \code{u74} - 1 if unemployed in 1974, 0 otherwise.
##'   \item \code{u75} - 1 if unemployed in 1975, 0 otherwise.
##'   \item \code{treat} - 1 if treated, 0 otherwise.
##' }
##'
##' @docType data
##' @keywords datasets
##' @name lalonde.psid
##' @usage data(lalonde.psid)
##' @format A data frame with 2675 rows and 12 variables
##' @references LaLonde, Robert J. (1986). Evaluating the Econometric
##' Evaluations of Training Programs with Experimental Data. The
##' American Economic Review, 76(4), 604--620.
NULL


##' Calculate sensitivity of causal estimates to unmeasured confounding.
##'
##' This function performs a sensitivity analysis of causal effects
##' different assumptions about unmeasured confounding, as described
##' by Blackwell (2013).
##'
##'
##' @param model.y outcome model object. Currently only handles
##' \code{lm} objects.
##' @param model.t propensity score model. Currently assumes a
##' \code{glm} object.
##' @param cov.form one-sided formula to describe any covariates to be
##' included in the parital R^2 calculations.
##' @param confound function that calculates the confounding
##' function. This function must take arguments \code{alpha},
##' \code{pscores}, and \code{treat}. Defaults to
##' \code{\link{one.sided}}. Other functions included with the package are
##' \code{\link{one.sided.att}}, \code{\link{alignment}}, and \code{\link{alignment.att}}.
##' @param data data frame to find the covariates from \code{cov.form}.
##' @param alpha vector of confounding values to pass the confounding
##' function. Defaults to 11 points from -0.5 to 0.5 for binary
##' outcome variable, and 11 points covering the
##' a interval with width equal to the inter-quartile range and centered at 0 for
##' non-binary outcome variables.
##' @return Returns an object of class \code{causalsens}.
##' \itemize{
##'   \item \code{sens} data frame containing alpha values, partial
##' R^2s, estimates, and 95% confidence intervals
##' \item \code{partial.r2} vector of partial R^2 values for the
##' covariates to compare to sensitivity analysis results.
##' }
##' @export
causalsens <- function(model.y, model.t, cov.form, confound = one.sided, data, alpha) {

  if (inherits(model.y, "glm")) {
    stop("Only works for linear outcome models right now. Check back soon.")
  }
  y.dat <- model.frame(model.y)
  t.dat <- model.frame(model.t)
  c.dat <- model.frame(cov.form, data)
  pscores <- fitted(model.t)
  rn.y <- row.names(y.dat)
  rn.t <- row.names(t.dat)
  t.name <- colnames(t.dat)[1]

  if (!identical(rn.y, rn.t)) {
    bothrows <- intersect(rn.y, rn.t)
    y.dat <- y.dat[bothrows,]
    t.dat <- t.dat[bothrows,]
    c.dat <- c.dat[bothrows,]
    pscores <- pscores[bothrows]
  }

  c.dat <- c.dat[,!(colnames(c.dat) %in% colnames(y.dat))]
  y.dat <- cbind(y.dat, c.dat)

  if (missing(alpha)) {
    if (length(unique(y.dat[,1])) == 2) {
      alpha <- seq(-0.5, 0.5, length = 11)
    } else {
      iqr <- quantile(y.dat[,1], 0.75) - quantile(y.dat[,1], 0.25)
      alpha <- seq(-iqr/2, iqr/2, length = 11)
    }
  }

  if ("(weights)" %in% colnames(y.dat)) {
    colnames(y.dat)[colnames(y.dat) == "(weights)"] <- as.character(model.y$call$weights)
  }

  rsq.form <- cov.form
  rsq.form[[2]] <- as.name("y.adj")
  rsq.form[[3]] <- cov.form[[2]]

  all.covs <- union(all.vars(model.y$terms[[3]]), all.vars(cov.form))
  all.covs <- all.covs[all.covs != t.name]
  partial.form <- as.formula(paste(colnames(y.dat)[1],
                                   paste(all.covs, collapse = " + "), sep = " ~ "))

  sens.form <- model.y$terms
  sens.form[[2]] <- as.name("y.adj")

  sens <- matrix(NA, nrow = length(alpha), ncol = 5)
  colnames(sens) <- c("rsqs","alpha", "estimate", "lower", "upper")
  sens[,"alpha"] <- alpha
  y.dat$y.adj <- NA
  for (j in 1:length(alpha)) {
    adj <- do.call(confound, list(alpha = alpha[j], pscores = pscores, treat = t.dat[,1]))
    y.dat$y.adj <- y.dat[,1] - adj
    s.out <- update(model.y, formula. = sens.form, data = y.dat)
    r.out <- update(model.y, formula = rsq.form, data = y.dat[t.dat[,1] == 0,])
    sens[j,4:5] <- confint(s.out)[t.name,]
    sens[j,3] <- coef(s.out)[t.name]
    sens[j,1] <- alpha[j]^2 * var(t.dat[,1])/var(residuals(s.out))
  }

  partial.out <- lm(partial.form, data = y.dat[t.dat[,1] == 0,])
  dropmat <- drop1(partial.out)
  prsqs <- dropmat[-1,2]/dropmat[-1,3]
  names(prsqs) <- rownames(dropmat)[-1]
  out <- list(sens = data.frame(sens), partial.r2 = prsqs)
  class(out) <- "causalsens"
  return(out)

}

##' @S3method print causalsens
print.causalsens <- function(x, ...) {
  cat("Sensitivity analysis output: \n")
  print(x$sens)
  cat("\nCovariate partial r-squared: \n")
  print(x$partial.r2)
  invisible()
}

##' @S3method summary causalsens
summary.causalsens <- function(object, ...) {
  cat("Maximum estimate: ", max(object$sens$estimate), "\n")
  cat("Minimum estimate: ", min(object$sens$estimate), "\n")
  cat("Parial R^2 range: [", min(object$sens$rsqs), ", ", max(object$sens$rsqs), "]\n", sep = "")
  cat("Covariate R^2 range: [", min(object$partial.r2), ", ", max(object$partial.r2), "]\n", sep = "")
}

##' Plot a causal sensitivity analysis.
##'
##' Plot the results of a sensitivity analysis against unmeasured
##' confounding as perfomed by \code{\link{causalsens}}
##' @S3method plot causalsens
##' @method plot causalsens
##' @export
##' @param x \code{causalsens} object.
##' @param type a string taking either the value \code{"r.squared"}
##' (default), which plots the estimated effects as a function of the
##' partial R-squared values, or \code{"raw"}, which plots them as a
##' function of the raw confounding values, \code{alpha}.
##' @param ... other parameters to pass to the plot.
plot.causalsens <- function(x, type = "r.squared", ...) {
  m <- match.call(expand.dots = TRUE)
  m[[1L]] <- quote(graphics::plot)
  if (m$type == "r.squared") {
    m$x <- sign(x$sens$alpha) * x$sens$rsqs
    if (is.null(m$xlab)) {
      m$xlab <- "Variance explained by confounding"
    }
  } else if (type == "raw") {
    m$x <- x$sens$alpha
    if (is.null(m$xlab)) {
      m$xlab <- "Amount of confounding"
    }
  } else {
    stop("type must be 'r.squared' or 'raw'")
  }
  if (is.null(m$ylim)) {
    m$ylim <- c(min(x$sens$lower), max(x$sens$upper))
  }
  if (is.null(m$ylab)) {
    m$ylab <- "Estimated effect"
  }
  m$y <- x$sens$estimate
  m$type <- "l"
  eval(m, parent.frame())
  ## plot(x = xpoints, y = x$sens$estimate, type = "l", ylim = ylim,
  ##      xlab = xlab, ylab = "Estimated effect", ...)
  abline(h = 0, col = "grey")
  abline(v = 0, col = "grey")
  polygon(x = c(m$x, rev(m$x)),
          y = c(x$sens$lower, rev(x$sens$upper)),
          col =rgb(0.5, 0.5, 0.5, alpha = 0.5), border = NA)
  if (type == "r.squared") {
    points(x = x$partial.r2, y = rep(0, length(x$partial.r2)), pch = 4)
    points(x = -x$partial.r2, y = rep(0, length(x$partial.r2)), pch = 4)
  }
  invisible()
}


##' Confounding functions
##'
##' Various confounding functions for use with
##' \code{\link{causalsens}}.
##' @aliases one.sided alignment one.sided.att alignment.att
##' @param alpha vector of confounding values to use in the
##' sensitivity analysis.
##' @param pscores vector of propensity scores for each unit.
##' @param treat vector of treatment values for each unit.
##' @export
one.sided <- function(alpha, pscores, treat) {
  adj <- alpha * (1 - pscores) * treat - alpha * pscores * (1 - treat)
  return(adj)
}
##' @rdname one.sided
##' @export
alignment <- function(alpha, pscores, treat) {
  adj <- alpha * (1 - pscores) * treat + alpha * pscores * (1 - treat)
  return(adj)
}
##' @rdname one.sided
##' @export
one.sided.att <- function(alpha, pscores, treat) {
  adj <- -alpha * pscores * (1 - treat)
  return(adj)
}
##' @rdname one.sided
##' @export
alignment.att <- function(alpha, pscores, treat) {
  adj <- alpha * pscores * (1 - treat)
  return(adj)
}
