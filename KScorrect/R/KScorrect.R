#' KScorrect: Lilliefors-Corrected Kolmogorov-Smirnoff Goodness-of-Fit Tests
#'
#' Implements the Lilliefors-corrected Kolmogorov-Smirnoff test for use in
#' goodness-of-fit tests.
#'
#' KScorrect implements the Lilliefors-corrected Kolmogorov-Smirnoff test for
#' use in goodness-of-fit tests, suitable when population parameters are unknown
#' and must be estimated by sample statistics. \emph{P}-values are estimated by
#' simulation. Coded to complement \code{\link[stats]{ks.test}}, it can be used
#' with a variety of continuous distributions, including normal, lognormal,
#' univariate mixtures of normals, uniform, loguniform, exponential, gamma, and
#' Weibull distributions.
#'
#' Functions to generate random numbers and calculate density, distribution, and
#' quantile functions are provided for use with the loguniform and mixture
#' distributions.
#'
#' @author Phil Novack-Gottshall \email{pnovack-gottshall@@ben.edu}
#' @author Steve C. Wang \email{scwang@@swarthmore.edu}
#' @name KScorrect-package
#' @aliases KScorrect-package KScorrect
#' @docType package
#'
#' @examples
#' # Get the package version and citation of KScorrect
#' packageVersion("KScorrect")
#' citation("KScorrect")
#'
#' x <- runif(200)
#' Lc <- LcKS(x, cdf="pnorm", nreps=999)
#' hist(Lc$D.sim)
#' abline(v = Lc$D.obs, lty = 2)
#' print(Lc, max=50)  # Print first 50 simulated statistics
#' # Approximate p-value (usually) << 0.05
#'
#' # Confirmation uncorrected version has increased Type II error rate when
#' #   using sample statistics to estimate parameters:
#' ks.test(x, "pnorm", mean(x), sd(x))   # p-value always larger, (usually) > 0.05
#'
#' x <- rlunif(200, min=exp(1), max=exp(10)) # random loguniform sample
#' Lc <- LcKS(x, cdf="plnorm")
#' Lc$p.value      # Approximate p-value: (usually) << 0.05
NULL
