#' @importFrom parallel mclapply
#' @importFrom assertthat assert_that is.string is.count are_equal
#' @importFrom R6 R6Class
NULL

#' Simulation based maximum likelihood estimation
#' 
#' @param model The model used for the estimation. 
#'   See \code{\link{create_jaatha_model}}.
#' @param data The data used for the estimation.
#'   See \code{\link{create_jaatha_data}}.
#' @param repetitions The number of independend optimizations that will be
#'   conducted. You should use a value greater than one here, to minimize
#'   the chance that the algorithms is stuck in a local maximum.
#' @param sim The number of simulations conducted for each step.
#' @param max_steps The maximal number of steps, in case Jaatha fails to 
#'   converge.
#' @param init_method Determines how the starting position of each repetition
#'   is chosen. See below for a description of the different options. 
#' @param cores The number of CPU cores that will be used for the simulations.
#'   The relies on the \pkg{parallel} package, and consequenlty only one
#'   core is supported on Windows.
#' @param sim_cache_limit The maximal number of simulations results that will be 
#'   cached. Cached results may be reused in following estimation steps if 
#'   they are within the current block. Reduce this value to save memory. 
#'   Setting this to a value smaller than \code{sim} disables caching. 
#' @param verbose If \code{TRUE}, information about the optimization algorithm
#'   is printed.
#' @param block_width The relative width of a block within jaatha will fit its
#'   local GLM. The default value is usually fine. Increasing this value may 
#'   help in case jaatha fails to converge, while you can try decreasing it if 
#'   the estimates of the likelihoods differ from the corrected values in the 
#'   'Correcting likelihoods for best estimates' phase.
#' @return A list contain the results. The list has the following entries:
#' \describe{
#'    \item{estimate}{The (approximated) maximum likelihood estimate}
#'    \item{loglikelihood}{The estimate log-likelihood of the estimate.}
#'    \item{converged}{A boolean indicating whether the optimization procedure
#'                     converged or not}
#'    \item{args}{The arguments provided to the jaatha function}
#' }
#' 
#' @section Initialization Methods:
#'   Jaatha has different options for determining the starting positions for 
#'   it's optimization procedure. The option \code{initial-search} will divide
#'   the parameter space in a number of equally sized block, estimate parameters
#'   within each block and use the estimates with the highest likelihood as
#'   starting positions. The option \code{zoom-in} starts with a block that
#'   is equal to the complete parameter space, estimate parameters in there,
#'   and then iteratively creates a smaller block around the estimates. Finally,
#'   \code{random} chooses random starting postions and
#'   \code{middle} will just start all repetitions at the middle of the 
#'   parameter space.
#'   
#' @author Paul Staab and Lisha Mathew
#' @export
jaatha <- function(model, data, 
                   repetitions = 3, 
                   sim = model$get_par_number() * 25, 
                   max_steps = 100, 
                   init_method = c("zoom-in", "initial-search", 
                                   "random", "middle"),
                   cores = 1,
                   verbose = TRUE,
                   sim_cache_limit = 10000,
                   block_width = 0.1) {
  
  # Check parameters
  assert_that(is_jaatha_model(model))
  assert_that(is_jaatha_data(data))
  assert_that(is.count(repetitions))
  assert_that(is.count(sim))
  assert_that(is.count(cores))
  assert_that(is.numeric(block_width) && length(block_width) == 1)
  assert_that(block_width > 0 && block_width < 1)
  
  # Setup
  log <- create_jaatha_log(model, data, repetitions, max_steps, verbose)
  if (sim_cache_limit < sim) sim_cache_limit <- 0
  sim_cache <- create_sim_cache(sim_cache_limit) #nolint
  
  # Get start positions
  log$log_initialization(init_method[1])
  start_pos <- get_start_pos(model, data, repetitions, sim, init_method, cores,
                             sim_cache = sim_cache, block_width = block_width)
  
  for (rep in 1:repetitions) {
    estimate <- start_pos[rep, ]
    log$log_new_rep(rep, estimate)
    likelihood <- -Inf
    last_lh_increase <- 0
    
    for (step in 1:max_steps) {
      block <- create_block(cbind(estimate - block_width * .5,
                                  estimate + block_width * .5), 
                            cut = TRUE)
      
      local_ml <- estimate_local_ml(block, model, data, sim, cores, sim_cache)
      if (is.null(local_ml)) {
        warning("The GLMs failed to converge. Aborting one repetition.")
        break
      }
      log$log_estimate(rep, step, local_ml)
      estimate <- local_ml$par
      
      if (local_ml$value > likelihood) {
        likelihood <- local_ml$value
        last_lh_increase <- step
      }
      
      if (step >= last_lh_increase + 10) {
        log$log_convergence(rep)
        break
      }
    }
  }
  
  # get presice llh values for best estimates
  log$log_llh_correction()
  best_values <- log$get_best_estimates(5)
  if (nrow(best_values) == 0) stop("No valid estimates.")
  for (i in 1:nrow(best_values)) {
    llh <- estimate_llh(model, data, as.numeric(best_values[i, -(1:3)]), #nolint 
                        100, cores, TRUE)
    log$log_estimate("final", i, llh, best_values[i, 3])
  }
  
  # return the results
  log$create_results()
}
