#' Run simulation for conservation schemes
#'
#' Run the metapopulation simulation for various conservation prioritization
#' schemes.
#'
#' @param weights A matrix of habitat weights. Each row corresponds to another
#'   scenario. Each column is a different habitat location.
#' @param reps Number of portfolios to simulate.
#' @param assess_freq The frequency (in generations) of spawner-recruit
#'   re-assessment. Passed to \code{\link{meta_sim}}.
#' @param burn Cycles to throw out as burn in.
#' @param risk_fn Type of variance or risk metric. By default takes the variance.
#'   Instead you can supply any function that takes a numeric vector and returns
#'   some single numeric value. E.g. CVaR.
#' @param ... Other values to pass to \code{\link{meta_sim}}.
#' @export
#' @return
#' Returns the portfolio mean and variance values and the simulation runs.

get_conserv_plans_mv <- function(weights, reps = 150, assess_freq = 5,
  burn = 1:30, risk_fn = var, ...) {
  n_pop = ncol(weights)
  port_mv <- list()
  port_out <- list()
  for(j in 1:nrow(weights)) {
    port_out[[j]] <- list()
    for(i in 1:reps) {
      port_out[[j]][[i]] <- meta_sim(b = weights[j, ], use_cache = FALSE,
        n_pop = n_pop, ...)
    }
    port_mv[[j]] <- plyr::ldply(port_out[[j]], function(x)
      get_port_vals(x, burn = burn, risk_fn = risk_fn))
  }
  return(list(port_mv = port_mv, port_out = port_out))
}
