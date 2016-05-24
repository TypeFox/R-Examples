#' Run a single metapopulation simulation.
#'
#' This is the master function for running \pkg{metafolio} simulations. It runs
#' a single iteration of a simulation. The arguments can be manipulated with
#' other functions in the package to use this function as part of a portfolio
#' analysis.
#'
#' @param n_t The number of years.
#' @param n_pop Number of populations
#' @param stray_decay_rate Rate that straying (exponentially) decays with
#'   distance.
#' @param stray_fraction Fraction of fish that stray from natal streams.
#' @param b Ricker density-dependent parameter. A vector with one numeric value
#'   per population.
#' @param spawners_0 A vector of spawner abundances at the start of the
#'   simulation. Length of the vector should equal the number of populations.
#' @param sigma_v Stock-recruit residual standard deviation of the
#'   log-deviations.
#' @param v_rho AR1 serial correlation of stock-recruit residuals.
#' @param a_width_param Width of the thermal curves by population.
#' @param optim_temp Optimal temperatures by population.
#' @param max_a Maximum Ricker productivity parameters (a) by population.  The
#'   value obtained at the optimum temperature. Note how the default argument
#'   uses the \code{\link{thermal_integration}} function.
#' @param env_type The type of environmental time series to generate. One of
#'   \code{"sine"}, \code{"arma"}, \code{"regime"}, \code{"linear"}, or
#'   \code{"constant"}. See \code{\link{generate_env_ts}}.
#' @param env_params Parameters to pass on to \code{\link{generate_env_ts}}.
#'   You must provide the appropriate list given your chosen type of
#'   environmental signal.
#' @param start_assessment Generation to start estimating the stock recruit
#'   relationship for escapement targets. The assessment is carried out using
#'   \code{\link{fit_ricker}}.
#' @param a_lim A vector of length two giving the lower and upper limits for
#'   Ricker a values. If a value is estimated beyond these limits it will be
#'   set to the limit value.
#' @param b_lim A vector of length two giving the lower and upper limits for
#'   the estimated Ricker b values *as fractions* of the previously assessed
#'   value. If a value is estimated beyond these limits it will be set to the
#'   limit value.
#' @param silence_warnings Should the warnings be skipped if the Ricker a or b
#'   values exceed their specified bounds? \code{meta_sim} will still print
#'   other warnings regardless of this argument value.
#' @param sigma_impl Implementation standard deviation for the implementation
#'   error beta distribution.
#' @param assess_freq How many generations before re-assessing Ricker a and b
#'   parameters.
#' @param use_cache Use the stochastically generated values (stock-recruit
#'   residuals and possibly environmental time series) from the previous run?
#'   See the Details section below.
#' @param cache_env Logical: Should the environmental time series be cached? If
#'   \code{use_cache = TRUE} then this will automatically happen. But, you
#'   could set \code{cache_env = TRUE} and \code{use_cache = FALSE} to only
#'   cache the environmental time series. See the Details section below.
#' @param add_straying Implement straying between populations?
#' @param add_impl_error Add implementation error? Implementation error is
#'   derived using \code{\link{impl_error}}.
#' @param decrease_b A numeric value to decrease all streams by each generation.
#'   This is intended to be used to simulate habitat loss, for example though
#'   stream flow reduction with climate change.
#' @param skip_saving_cache Logical: if \code{TRUE} then no data will be cached
#'   for the next iteration. This will save time when running many simulations.
#' @param debug Logical: if \code{TRUE} then \code{meta_sim} will print a number
#'   of debugging statements while it runs.
#' @details
#' To use either of the caching options, you must have run \code{meta_sim} at
#' least once in the current session with both caching arguments set to
#' \code{FALSE} to generate the cached values first. If you're running many
#' iterations of \code{meta_sim} and you want to cache, then the first iteration
#' should have both cache arguments set to \code{FALSE}, and subsequent runs can
#' set one or both to \code{TRUE}. Internally, \code{meta_sim} caches by writing
#' the appropriate data to an \code{.rda} file in a temporary directory.
#' @return
#' A list is returned that contains the following elements. All matrices that
#' are returned (except the straying matrix) feature populations along the
#' columns and generations/years along the rows.
#' \describe{
#'    \item{\code{A}}{A matrix of abundances.}
#'    \item{\code{F}}{A matrix of fishing mortality in numbers.}
#'    \item{\code{E}}{A matrix of realized escapement.}
#'    \item{\code{Eps}}{A matrix of (log) spawner-return residuals. These have been
#'    log-normal bias corrected so their expected value after exponentiation
#'    will be one.}
#'    \item{\code{A_params}}{A matrix of actual Ricker a parameters.}
#'    \item{\code{Strays_leaving}}{A matrix of strays leaving.}
#'    \item{\code{Strays_joining}}{A matrix of strays joining.}
#'    \item{\code{env_ts}}{A vector of the environmental time series.}
#'    \item{\code{stray_mat}}{The straying matrix. These fractions are constant across
#'    generations/years. Rows and columns are populations.}
#'    \item{\code{n_pop}}{The total possible populations as input in the simulation.}
#'    \item{\code{n_t}}{The number of generations/years the simulation was run for.}
#'    \item{\code{b}}{The original Ricker b values as specified.}
#'    \item{\code{Est_a}}{A matrix of estimated Ricker a values.}
#'    \item{\code{Est_b}}{A matrix of estimated Ricker b values.}
#' }
#' @export
#' @examples
#' arma_env_params <- list(mean_value = 16, ar = 0.1, sigma_env = 2, ma = 0)
#' base1 <- meta_sim(n_pop = 10, env_params = arma_env_params,
#'   env_type = "arma", assess_freq = 5)
#'
#' plot_sim_ts(base1, years_to_show = 70, burn = 30)

meta_sim <- function(
  n_t = 130,
  n_pop = 10,
  stray_decay_rate = 0.1,
  stray_fraction = 0.02,
  b = rep(1000, n_pop),
  spawners_0 = round(b),
  sigma_v = 0.7,
  v_rho = 0.4,
  a_width_param = c(seq(0.08, 0.04, length.out = n_pop/2), rev(seq(0.08, 0.04,
        length.out = n_pop/2))),
  optim_temp = seq(13, 19, length.out = n_pop),
  max_a = thermal_integration(n_pop),
  env_type = c("sine", "arma", "regime", "linear", "constant"),
  env_params = list(amplitude = 3.2, ang_frequency = 0.2, phase = runif(1,
      -pi, pi), mean_value = 15, slope = 0, sigma_env = 0.30),
  start_assessment = 20,
  a_lim = c(0.02, 4),
  b_lim = c(0.5, 1.5),
  silence_warnings = TRUE,
  sigma_impl = 0.1,
  assess_freq = 10,
  use_cache = FALSE,
  cache_env = FALSE,
  add_straying = TRUE,
  add_impl_error = TRUE,
  skip_saving_cache = FALSE,
  decrease_b = 0,
  debug = FALSE
  ) {

  if (use_cache | cache_env) {
    load(paste0(tempdir(), "/env_ts.rda"))
  } else {
    env_type <- env_type[1]
    env_ts <- switch(env_type,
      sine = generate_env_ts(n_t = n_t, type = "sine", sine_params =
        env_params),
      arma = generate_env_ts(n_t = n_t, type = "arma", arma_params =
        env_params),
      regime = generate_env_ts(n_t = n_t, type = "regime",
        regime_params = env_params),
      linear = generate_env_ts(n_t = n_t, type = "linear",
        linear_params = env_params),
      linear_arma = generate_env_ts(n_t = n_t, type = "linear_arma",
        linear_arma_params = env_params),
      constant = generate_env_ts(n_t = n_t, type = "constant",
        constant_params = env_params)
      )
    if(!skip_saving_cache) {
      save(env_ts, file = paste0(tempdir(), "/env_ts.rda"))
    }
  }

  # create vectors and matrices that are not developed iteratively:
  if (use_cache) {
    load(paste0(tempdir(), "/sim_dat.rda"))
  } else {

    # Figure out alpha parameters before running through the loops:
    A_params <- matrix(ncol = n_pop, nrow = n_t) # a parameters from Ricker
    for(j in 1:n_pop) {
      A_params[, j] <- thermal_curve_a(env_ts, optim_temp =
        optim_temp[j], max_a = max_a[j], width_param =
        a_width_param[j])
    }
    A_params[A_params<0] <- 0

    stray_mat <- generate_straying_matrix(n_pop = n_pop, stray_fraction
      = stray_fraction, stray_decay_rate = stray_decay_rate)

    epsilon_mat <- matrix(ncol = n_pop, nrow = n_t)
    epsilon_mat[1, ] <- rnorm(n_pop, mean = 0, sd = sigma_v)
    for(i in 2:n_t) {
      epsilon_mat[i, ] <- epsilon_mat[i - 1, ] * v_rho + rnorm(n_pop, mean =
        0 - (sigma_v^2)/2, sd = sigma_v) # stock-recruit residuals; note bias correction
    }
    # now develop random escapement targets at start of open access
    r_escp_goals <- matrix(nrow = start_assessment, ncol = n_pop, data =
      runif(n_pop*start_assessment, 0.1, 0.9))
    if(!skip_saving_cache) {
      save(stray_mat, epsilon_mat, r_escp_goals, A_params, file = paste0(tempdir(), "/sim_dat.rda"))
    }
  }

  assess_years <- seq(start_assessment, n_t, assess_freq)

  out <- metasim_base(n_pop = n_pop, n_t = n_t, spawners_0 = spawners_0, b = b,
    epsilon_mat = epsilon_mat, A_params = A_params, add_straying = add_straying,
    stray_mat = stray_mat, assess_years = assess_years, r_escp_goals =
    r_escp_goals, sigma_impl = sigma_impl, add_impl_error = add_impl_error,
    decrease_b = decrease_b, debug = debug)
  out$env_ts <- env_ts
  return(out)
}
