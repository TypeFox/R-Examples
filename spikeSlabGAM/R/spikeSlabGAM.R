#' @include ssGAMDesign.R
#' @include spikeAndSlab.R
{}

#' Generate posterior samples for a GAMM with spike-and-slab priors
#'
#' This function fits structured additive regression models with spike-and-slab
#' priors via MCMC. The spike-and-slab prior results in an SSVS-like term
#' selection that can be used to aid model choice, e.g. to decide between linear
#' and smooth modelling of numeric covariates or which interaction effects are
#' relevant. Model terms can be factors (\code{\link{fct}}), linear/polynomial
#' terms (\code{\link{lin}}),  uni- or bivariate splines (\code{\link{sm}},
#' \code{\link{srf}}), random intercepts (\code{\link{rnd}}) or Markov random
#' fields (\code{\link{mrf}}) and their interactions, i.e. an interaction
#' between two smooth terms yields an effect surface, an interaction between a
#' linear term and a random intercept yields random slopes, an interaction
#' between a linear term and a smooth term yields a varying coefficient term
#' etc.
#'
#' Implemented types of terms: \describe{ \item{\code{\link{u}}}{unpenalized
#' model terms associated with a flat prior (no selection)}
#' \item{\code{\link{lin}}}{linear/polynomial terms}
#' \item{\code{\link{fct}}}{factors} \item{\code{\link{sm}}}{univariate smooth
#' functions} \item{\code{\link{rnd}}}{random intercepts}
#' \item{\code{\link{mrf}}}{Markov random fields}
#' \item{\code{\link{srf}}}{surface estimation} } Terms in the formula that are
#' not in the list of specials (i.e. the list of term types above) are
#' automatically assigned an appropriate special wrapper, i.e. numerical
#' covariates \code{x} are treated as \code{\link{lin}(x) + \link{sm}(x)},
#' factors \code{f} (and numerical covariates with very few distinct values, see
#' \code{\link{ssGAMDesign}}) are treated as \code{\link{fct}(f)}. Valid model
#' formulas have to satisfy the following conditions: \enumerate{ \item All
#' variables that are involved in interactions have to be present as main
#' effects as well, i.e. {\code{y ~ x1 + x1:x2}} will produce an error because
#' the main effect of \code{x2} is missing. \item Interactions between
#' unpenalized terms not under selection (i.e. terms of type \code{\link{u}})
#' and penalized terms are not allowed, i.e. \code{y ~ u(x1)* x2} will produce an
#' error. \item If main effects are specified as special terms, the interactions
#' involving them have to be given as special terms as well, i.e. \code{y ~
#' lin(x1) + lin(x2) + x1:x2} will produce an error. }
#'
#' The default prior settings for Gaussian data work best for a response with
#' unit variance. If your data is scaled very differently, either rescale the
#' response (recommended) or adjust the hyperparameters accordingly.
#'
#' Sampling of the chains is done in parallel using package \code{multicore} if
#' available. If not, a socket cluster set up with \code{snow} is used where
#' available. Use \code{options(mc.cores = foo)} to set the (maximal) number of
#' processes spawned by the parallelization. If \code{options()$mc.cores} is
#' unspecified, snow uses 8.
#'
#' @param formula the model formula (see details below and
#'   \code{\link{ssGAMDesign}}).
#' @param data the data containing the variables in the function
#' @param ... additional arguments for \code{\link{ssGAMDesign}}
#' @param family (character) the family of the response, defaults to
#'   normal/Gaussian response, \code{"poisson"} and \code{"binomial"} are
#'   implemented as well.
#' @param hyperparameters A list of arguments specifying the prior settings. See
#'   \code{\link{spikeAndSlab}}.
#' @param model A list of arguments describing the model structure. See
#'   \code{\link{spikeAndSlab}}. User-supplied \code{groupIndicators} and
#'   \code{H} entries will be overwritten by \code{\link{ssGAMDesign}}.
#' @param mcmc A list of arguments specifying MCMC sampler options. See
#'   \code{\link{spikeAndSlab}}.
#' @param start A list of starting values for the MCMC sampler. See
#'   \code{\link{spikeAndSlab}}. Use \code{start = list(seed =<YOUR_SEED>)} to set
#'   the RNG seed for reproducible results.
#' @return an object of class \code{spikeSlabGAM} with methods
#'   \code{\link{summary.spikeSlabGAM}}, \code{\link{predict.spikeSlabGAM}}, and
#'   \code{\link{plot.spikeSlabGAM}}.
#' @seealso \code{\link{ssGAMDesign}} for details on model specification,
#'   \code{\link{spikeAndSlab}} for more details on the MCMC sampler and prior
#'   specification, and \code{\link{ssGAM2Bugs}} for MCMC diagnostics. Check out
#'   the vignette for theoretical background and code examples.
#' @author Fabian Scheipl
#' @references Fabian Scheipl (2011). \code{spikeSlabGAM}: Bayesian Variable
#'   Selection, Model Choice and Regularization for Generalized Additive Mixed
#'   Models in R. \emph{Journal of Statistical Software}, \bold{43}(14), 1--24.
#'
#'   Fabian Scheipl, Ludwig Fahrmeir, Thomas Kneib (2012). Spike-and-Slab Priors
#'   for Function Selection in Structured Additive Regression Models.
#'   \emph{Journal of the American Statistical Association}, \bold{107}(500),
#'   1518--1532.
#' @import utils
#' @export
#' @examples
#' \dontrun{
#' ## examples not run due to time constraints on CRAN checks.
#' ## full examples below should take about 2-4 minutes.
#'
#' set.seed(91179)
#' n <- 400
#' d <- data.frame(f1 = sample(gl(3, n/3)), f2 = sample(gl(4,
#' 						n/4)), x1 = runif(n), x2 = runif(n), x3 = runif(n))
#' # true model:
#' #   - interaction between f1 and x1
#' #   - smooth interaction between x1 and x2,
#' #   - x3 and f2 are noise variables without influence on y
#' nf1 <- as.numeric(d$f1)
#' d$f <- with(d, 5 * (nf1 + 2 * sin(4 * pi * (x1 - 0.2) *
#' 									(x2 - 0.7)) - nf1 * x1))
#' d$y <- with(d, scale(f + rnorm(n)))
#' d$yp <- with(d, rpois(n, exp(f/10)))
#'
#' # fit & display the model:
#' m1 <- spikeSlabGAM(y ~ x1 * f1 + f1 * f2 + x3 * f1 +
#' 				x1 * x2, data = d)
#' summary(m1)
#'
#' # visualize estimates:
#' plot(m1)
#' plot(m1, cumulative = FALSE)
#' (plotTerm("fct(f1):fct(f2)", m1, aggregate = median))
#' print(p <- plotTerm("sm(x1):sm(x2)", m1, quantiles = c(0.25,
#' 						0.75), cumulative = FALSE))
#'
#' # change MCMC settings and priors:
#' mcmc <- list(nChains = 3, burnin = 100, chainLength = 1000,
#' 		thin = 3, reduceRet = TRUE)
#' hyper <- list(gamma = c(v0 = 5e-04), tau2 = c(10,
#' 				30), w = c(2, 1))
#'
#' # complicated formula example, poisson response:
#' m2 <- spikeSlabGAM(yp ~ x1 * (x2 + f1) + (x2 + x3 + f2)^2 -
#'          sm(x2):sm(x3), data = d,
#'   				family = "poisson", mcmc = mcmc,
#' 		      hyperparameters = hyper)
#' summary(m2)
#' plot(m2)
#'
#' # quick&dirty convergence diagnostics:
#' print(b <- ssGAM2Bugs(m2))
#' plot(b)
#' }
spikeSlabGAM <- function(formula,
  data,
  ...,
  family ="gaussian",
  hyperparameters = list(),    # prior parameters
  model = list(),              # model structure
  mcmc = list(),               # Gibbs sampling options
  start = list()               # start values for the Gibbs sampler
)
{
  #require(RegularBayes)
  D <- ssGAMDesign(formula = formula, data = data, ...)
  model$groupIndicators <- D$groupIndicators
  model$H <- D$H
  res <- spikeAndSlab(family = family,
    y = D$response,           # response (n x 1)
    X = D$Design,             # design matrix with covariates (n x q)
    hyperparameters = hyperparameters,    # prior parameters
    model = model,              # model structure
    mcmc = mcmc,      	      # Gibbs sampling options
    start = start)              # start values for the Gibbs sampler
  res$formula <- D$formula
  res$data <- data
  res$predvars <- D$predvars
  class(res) <- c("spikeSlabGAM")
  return(res)
}



