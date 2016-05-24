#' Monte Carlo asset weights into portfolios
#'
#' Monte Carlo the asset weights into portfolios and record the simulation
#' output and portfolio metrics (mean and variance).
#'
#' @param weights_matrix A matrix of asset weights. The columns correspond to
#'   the different assets and the rows correspond to the simulation iterations.
#' @param n_sims The number of simulations to run.
#' @param mean_b The mean Ricker capacity value.
#' @param burn The number of years to discard as burn in.
#' @param ... Anything else to pass to \code{\link{meta_sim}}.
#' @return A list object with three elements: \code{port_vals} (a matrix with a
#'   column of mean rate of change and variance of rate of change),
#'   \code{n_sims} (the number of simulations ran), and \code{sims_out} (a list
#'   in which each element corresponds to the output from the run of
#'   \code{\link{meta_sim}}.
#' @seealso \code{\link{meta_sim}}, \code{\link{create_asset_weights}}
#' @export
#' @examples
#' weights_matrix <- create_asset_weights(n_pop = 4, n_sims = 3,
#'   weight_lower_limit = 0.001)
#' mc_ports <- monte_carlo_portfolios(weights_matrix = weights_matrix,
#'   n_sims = 3, mean_b = 1000)

monte_carlo_portfolios <- function(weights_matrix, n_sims = 500,
  mean_b = 1000, burn = 1:30, ...) {

  b_vals_matrix <- weights_matrix * mean_b * ncol(weights_matrix)
  n_pop <- ncol(weights_matrix)
  port_vals = matrix(ncol = 2, nrow = n_sims)
  sims_out <- list() # to store simulation output

  for(k in 1:n_sims) {
    b_vals <- b_vals_matrix[k, ]
    if(k == 1)
      sims_out[[k]] <- meta_sim(b = b_vals, use_cache = FALSE, n_pop =
        n_pop, ...)
    else
      sims_out[[k]] <- meta_sim(b = b_vals, use_cache = TRUE, n_pop =
        n_pop, ...)
    port.x <- rowSums(sims_out[[k]]$A[-burn, ])
    ret.x <- diff(log(port.x))
    var.x <- var(ret.x)
    mean.x <- mean(ret.x)
    port_vals[k, ] <- c(mean.x, var.x)
    if(k%%20 == 0 | k == n_sims) print(paste0("Done ", k, " of ", n_sims, "."))
  }
  return(list(port_vals = port_vals, n_sims = n_sims, sims_out = sims_out))
}
