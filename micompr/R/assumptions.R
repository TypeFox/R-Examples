# Copyright (c) 2016 Nuno Fachada
# Distributed under the MIT License (http://opensource.org/licenses/MIT)

#' Parametric tests assumptions
#'
#' Generic function to get the assumptions for parametric tests applied to the
#' comparison of output observations.
#'
#' @param obj Object from which to get the assumptions.
#'
#' @return Assumptions for parametric tests applied to the comparison of
#' outputs.
#'
#' @export
#'
#' @seealso \code{\link{assumptions.cmpoutput}},
#' \code{\link{assumptions.micomp}}
#'
assumptions <- function(obj) UseMethod("assumptions")

#' Determine the assumptions for the MANOVA test
#'
#' Determine two assumptions for the MANOVA test: a) multivariate normality
#' of each group; b) homogeneity of covariance matrices.
#'
#' @param data Data used for the MANOVA test (rows correspond to observations,
#' columns to dependent variables).
#' @param factors Groups to which rows of \code{data} belong to (independent
#' variables).
#'
#' @return An object of class \code{assumptions_manova} which is a list
#' containing two elements:
#' \describe{
#'  \item{\code{mvntest}}{List of results from the Royston multivariate
#'        normality test (\code{\link[MVN]{roystonTest}}), one result per
#'        group.}
#'  \item{\code{vartest}}{Result of Box's M test for homogeneity of covariance
#'        matrices (\code{\link[biotools]{boxM}}).}
#' }
#'
#' @export
#'
#' @note This function requires the \code{MVN} and \code{biotools} packages.
#'
#' @examples
#'
#' # Determine the assumptions of applying MANOVA to the iris data
#' # (i.e. multivariate normality of each group and homogeneity of covariance
#' # matrices)
#' a <- assumptions_manova(iris[, 1:4], iris[, 5])
#'
assumptions_manova <- function(data, factors) {

  # Don't return MANOVA assumptions if the required packages are not installed.
  if (!requireNamespace("MVN", quietly = TRUE) ||
      !requireNamespace("biotools", quietly = TRUE) ) {
    message("MANOVA assumptions require 'MVN' and 'biotools' packages.")
    return(NULL)
  }

  # `assumpt` will be a list containing the test for multivariate normality
  # and the test for the homogeneity of covariance matrices
  assumpt <- list()
  class(assumpt) <- "assumptions_manova"
  assumpt$mvntest <- list()
  assumpt$vartest <- NULL

  # Number of variables
  nvars <- dim(data)[2]

  # Cycle through factors for the multivariate normality tests
  for (f in unique(factors)) {

    # Indexes of current factor level
    fidx <- factors == f
    nobs <- sum(fidx)

    # Royston test requires at least four observations
    if (nobs > 3) {

      # Royston test requires that there are more observations than variables
      # for each group
      if (nobs > nvars) {

        assumpt$mvntest[[f]] <-
          MVN::roystonTest(data[fidx, ])

      } else {

        # If there are no more observations than variables for current group,
        # then perform test with less variables and warn the user
        assumpt$mvntest[[f]] <-
          MVN::roystonTest(data[fidx, 1:min(nobs - 1, nvars)])
        warning(paste("Royston test requires more observations than ",
                      "(dependent) variables (DVs). Reducing number of ",
                      "variables from ", nvars, " to ", nobs - 1," in group '",
                      f, "'.", sep = ""), call. = F)
      }
    } else {

      # Don't perform test if number of observations is less than 4
      assumpt$mvntest[[f]] <- NA
      warning(paste("Royston test requires at least 4 observations ",
                    "(independent variables), but there are only ", nobs,
                    " observations in group '", f, "'. Test not performed.",
                    sep = ""), call. = F)
    }
  }

  # Perform the homogeneity of covariance matrices test (Box's M)
  maxvars <- min(table(factors) - 1, nvars)
  assumpt$vartest <- biotools::boxM(data[, 1:maxvars], factors)

  assumpt
}

#' Determine the assumptions for the parametric comparison test
#'
#' Determine two assumptions for the parametric comparison tests (i.e. either
#' \code{\link[stats]{t.test}} or \code{\link[stats]{aov}}) for each principal
#' component, namely: a) univariate normality of each group; b) homogeneity of
#' variances.
#'
#' @param data Data used in the parametric test (rows correspond to
#' observations, columns to principal components).
#' @param factors Groups to which rows of \code{data} belong to.
#'
#' @return An object of class \code{assumptions_paruv} which is a list
#' containing two elements:
#' \describe{
#'  \item{\code{uvntest}}{List of results from the Shapiro-Wilk normality test
#'        (\code{\link[stats]{shapiro.test}}), one result per group per
#'        principal component.}
#'  \item{\code{vartest}}{Result of Bartlett test for homogeneity of variances
#'        (\code{\link[stats]{bartlett.test}}).}
#' }
#'
#' @export
#'
#' @examples
#'
#' # Determine the assumptions of applying ANOVA to each column (dependent
#' # variable) of the iris data (i.e. normality of each group and homogeneity of
#' # variances)
#' a <- assumptions_paruv(iris[, 1:4], iris[, 5])
#'
assumptions_paruv <- function(data, factors) {

  # Number of variables
  nvars <- dim(data)[2]
  if (is.null(nvars)) nvars <- 1

  # `assumpt` will be a list containing the tests for univariate normality
  # assumptions and the test for the homogeneity of variances assumption.
  assumpt <- list()
  class(assumpt) <- "assumptions_paruv"
  assumpt$uvntest <- list()
  assumpt$vartest <- list()

  # Cycle through each variable
  for (d in 1:nvars) {

    # Get current data
    if (nvars == 1) {
      currdata <- data
    } else {
      currdata <- data[, d]
    }

    # Perform univariate normality tests for each group for current variable
    for (f in unique(factors)) {
      if (!is.list(assumpt$uvntest[[f]])) assumpt$uvntest[[f]] <- list()
      assumpt$uvntest[[f]][[d]] <- stats::shapiro.test(currdata[factors == f])
    }

    # Perform the homogeneity of variances test for current variable
    assumpt$vartest[[d]] <- stats::bartlett.test(currdata ~ factors)

  }

  assumpt

}

#' Print information about the assumptions of the MANOVA test
#'
#' Print information about objects of class \code{assumptions_manova}, which
#' represent the assumptions of the MANOVA test performed on a comparison of
#' outputs.
#'
#' @param x Object of class \code{assumptions_manova}.
#' @param ... Currently ignored.
#'
#' @return The argument \code{x}, invisibly, as for all \code{\link{print}}
#' methods.
#'
#' @export
#'
#' @examples
#'
#' # Print information concerning the assumptions of applying MANOVA to the iris
#' # data (i.e. multivariate normality of each group and homogeneity of
#' # covariance matrices)
#' assumptions_manova(iris[, 1:4], iris[, 5])
#'
print.assumptions_manova <- function(x, ...) {

  cat("Royston test (Multivariate Normality):\n")
  for (grp in names(x$mvntest)) {
    if (methods::is(x$mvntest[[grp]], "royston")) {
      cat("\tP-value '", grp, "': ",
          x$mvntest[[grp]]@p.value, "\n", sep = "")
    } else {
      cat("\tTest not performed.\n")
    }
  }
  cat("Box's M test (Homogeneity of Covariance Matrices):\n")
  cat("\tP-value:", x$vartest$p.value, "\n")

  invisible(x)
}

#' Print information about the assumptions of the parametric test
#'
#' Print information about objects of class \code{assumptions_paruv}, which
#' represent the assumptions of the parametric test (i.e. either
#' \code{\link[stats]{t.test}} or \code{\link[stats]{aov}}) performed on a
#' comparison of outputs.
#'
#' @param x Object of class \code{assumptions_paruv}.
#' @param ... Currently ignored.
#'
#' @return The argument \code{x}, invisibly, as for all \code{\link{print}}
#' methods.
#'
#' @export
#'
#' @examples
#'
#' # Print information about the assumptions of applying ANOVA to each column
#' # (dependent variable) of the iris data (i.e. normality of each group and
#' # homogeneity of variances)
#' assumptions_paruv(iris[, 1:4], iris[, 5])
#'
print.assumptions_paruv <- function(x, ...) {

  maxvars <- min(5, length(x$uvntest[[1]]))

  cat("Shapiro-Wilk test (Normality):\n")
  for (grp in names(x$uvntest)) {
    cat("\tP-value(s) '", grp, "': ", sep = "")
    for (i in 1:maxvars) {
      cat("", x$uvntest[[grp]][[i]]$p.value)
    }
    if (length(x$uvntest[[1]]) > 5) {
      cat(" ... \n")
    } else {
      cat("\n")
    }
  }

  cat("Bartlett test (Homogeneity of Variances):\n\tP-value(s): ")
  for (i in 1:maxvars) {
    cat("", x$vartest[[i]]$p.value)
  }
  if (length(x$uvntest[[1]]) > 5) {
    cat(" ... \n")
  } else {
    cat("\n")
  }

  invisible(x)
}

#' Plot \emph{p}-values for testing the multivariate normality assumptions of
#' the MANOVA test
#'
#' Plot method for objects of class \code{\link{assumptions_manova}} which
#' presents a bar plot containing the \emph{p}-values produced by the Royston
#' multivariate normality test (\code{\link[MVN]{roystonTest}}) for each group
#' being compared.
#'
#' @param x Objects of class \code{\link{assumptions_manova}}.
#' @param ... Extra options passed to \code{\link[graphics]{barplot}}. The
#' \code{col} parameter defines colors for \emph{p}-values below 1, 0.05 and
#' 0.01, respectively.
#'
#' @return None.
#'
#' @export
#'
#' @importFrom graphics plot
#'
#' @examples
#'
#' # Plot the Royston test p-value for multivariate normality of each group
#' # (species) of the iris data
#' plot(assumptions_manova(iris[, 1:4], iris[, 5]))
#'
#' # Plot the same data with logarithmic scale for p-values
#' plot(assumptions_manova(iris[, 1:4], iris[, 5]), log = "y")
#'
plot.assumptions_manova <- function(x, ...) {

  # Was a color specified?
  params <- list(...)
  if (exists("col", where = params)) {
    col <- params$col
    params$col <- NULL
  } else {
    col <- c("darkgreen", "yellow", "red")
  }

  # Get the p-values to plot
  pvals <- sapply(x$mvntest, function(x) x@p.value)

  # Plot the p-values in a bar plot
  params$height <- pvals
  params$main <- sprintf("Royston test p-values (%d PCs)",
                         dim(x$mvntest[[1]]@dataframe)[2])
  params$sub <- "Multivariate normality"
  params$xlab <- "Groups"
  params$ylab <- "Probability"
  params$col <- pvalcol(pvals, col)
  do.call("barplot", params)

  invisible(NULL)
}

#' Plot \emph{p}-values for testing the assumptions of the parametric tests used
#' in output comparison
#'
#' Plot method for objects of class \code{\link{assumptions_paruv}} containing
#' \emph{p}-values produced by testing the assumptions of the parametric tests
#' used for comparing outputs.
#'
#' One bar plot is presented for the Bartlett test
#' (\code{\link[stats]{bartlett.test}}), showing the respective \emph{p}-values
#' along principal component. \emph{s} bar plots are presented for the
#' Shapiro-Wilk (\code{\link[stats]{shapiro.test}}), where \emph{s} is the
#' number of groups being compared; individual bars in each plot represent the
#' \emph{p}-values associated with each principal component.
#'
#' @param x Objects of class \code{\link{assumptions_paruv}}.
#' @param ... Extra options passed to \code{\link[graphics]{barplot}}. The
#' \code{col} parameter defines colors for \emph{p}-values below 1, 0.05 and
#' 0.01, respectively.
#'
#' @return None.
#'
#' @export
#'
#' @importFrom graphics plot
#'
#' @examples
#'
#' # Plot the Shapiro-Wilk and Bartlett test p-values for each dependent
#' # variable of the iris data
#' plot(assumptions_paruv(iris[, 1:4], iris[, 5]))
#'
#' # Plot the same data with logarithmic scale for p-values
#' plot(assumptions_paruv(iris[, 1:4], iris[, 5]), log = "y")
#'
plot.assumptions_paruv <- function(x, ...) {

  # Was a color specified?
  params <- list(...)
  if (exists("col", where = params)) {
    col <- params$col
    params$col <- NULL
  } else {
    col <- c("darkgreen", "yellow", "red")
  }

  # Number of vars in the PC plots
  nvars <- length(x$uvntest[[1]])

  # Plot the Bartlett test p-values by PC
  vardata <- sapply(x$vartest, function(x) x$p.value)
  params_bart <- params
  params_bart$height <- vardata
  params_bart$names.arg <- as.character(1:nvars)
  params_bart$main <- "p-values for the Bartlett test"
  params_bart$sub <- "Homogeneity of Variances"
  params_bart$xlab <- "PC"
  params_bart$ylab <- "Probability"
  params_bart$col <- pvalcol(vardata, col)
  do.call("barplot", params_bart)

  # Plot the Shapiro-Wilk p-values by PC for each factor
  for (grp in names(x$uvntest)) {
    normdata <- sapply(x$uvntest[[grp]], function(x) x$p.value)
    params_sw <- params
    params_sw$height <- normdata
    params_sw$names.arg <- as.character(1:nvars)
    params_sw$sub <- grp
    params_sw$main <- "p-values for the SW normality test"
    params_sw$xlab <- "PC"
    params_sw$ylab <- "Probability"
    params_sw$col <- pvalcol(normdata, col)
    do.call("barplot", params_sw)
  }

  invisible(NULL)

}
