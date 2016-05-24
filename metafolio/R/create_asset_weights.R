#' Create an asset weights matrix
#'
#' Create an asset weight matrix to run through the Monte Carlo
#' algorithm and test possible portfolios.
#'
#' @param n_pop The number of subpopulations.
#' @param n_sims The number of simulations.
#' @param weight_lower_limit The lowest fraction allowed for a subpopulation
#'   weight. For example, a value of 0.02 means a subpopulation will at least
#'   be assigned 2\% of the total capacity
#' @examples
#' create_asset_weights(n_pop = 5, n_sims = 10, weight_lower_limit = 0.001)
#' @export
#' @return A matrix. The columns represent subpopulations. The rows
#' represent simulation repetitions.
create_asset_weights <- function(n_pop, n_sims,
  weight_lower_limit = 0.02) {
  weights_matrix <- matrix(ncol = n_pop, nrow = n_sims)
  for (i in 1:n_sims) {
    w_i <- 0
    while (min(w_i) < weight_lower_limit) { # ensure no weights are too small
      w_i <- runif(n_pop)
      w_i <- w_i / sum(w_i) # weights that add to one
    }
    weights_matrix[i,] <-  w_i
  }
  weights_matrix
}

