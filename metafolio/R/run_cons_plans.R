#' Run conservation plans and return the portfolio mean and variance values
#'
#' This function takes a set of weights representing different conservation
#' plans and gets the mean and variance in portfolio space. This function allows
#' a maximally complicated set of weights to accommodate all possible scenarios.
#' It can accommodate different spatial strategies of conservation, conserving
#' different numbers of populations, and a lack of knowledge. You can do this by
#' how you set your \code{w} weight object. See the example.
#'
#' @param w A (nested) list of weights. The first list level contains the
#'   different plans. The next level contains repetitions for a given plan.
#'   E.g. \code{cp[[2]][[1]]} contains the first iteration of the second
#'   conservation plan. Each end element should be a matrix of weights with one
#'   row and the number of columns equal to the number of subpopulations.
#' @param env_type The environmental type to pass to
#'   \code{\link{generate_env_ts}}
#' @param env_params The environmental parameters to pass to
#'   \code{\link{generate_env_ts}}
#' @param show_progress Logical: show an indication of progress?
#' @param burn Cycles to throw out as burn in
#' @param assess_freq How frequently (in years) to re-assess the Ricker a and b
#'   values.
#'  @param risk_fn A risk function to use. Can be any function that takes a
#'    numeric vector and returns a single value. Suggested values include
#'    \code{var}, or \code{\link{VaR}}, or \code{\link{CVaR}}. Defaults to
#'    variance.
#' @param ... Other values to pass to \code{\link{meta_sim}}
#' @export
#' @return A list with two high-level elements: the mean variance output
#'   (\code{plans_mv}) and the raw simulation output (\code{plans_port}).
#'   Within \code{plans_mv}, each element of the list contains a conservation
#'   plan. Each row of the data frames represents a trial run. Within
#'   \code{plans_port}, each first level of the list contains a weight element
#'   and each second level of the list contains a replicate.
#' @examples
#' \dontrun{
#' set.seed(1)
#' w_plans <- list()
#' w_plans[[1]] <- c(5, 1000, 5, 1000, 5, 5, 1000, 5, 1000, 5)
#' w_plans[[2]] <- c(5, 5, 5, 1000, 1000, 1000, 1000, 5, 5, 5)
#' w_plans[[3]] <- c(rep(1000, 4), rep(5, 6))
#' w_plans[[4]] <- rev(w_plans[[3]])
#' plans_name_sp <- c("Full range of responses", "Most stable only",
#' "Lower half", "Upper half")
#'  n_trials <- 50 # number of trials at each n conservation plan
#'  n_plans <- 4 # number of plans
#'  num_pops <- c(2, 4, 8, 16) # n pops to conserve
#'  w <- list()
#'  for(i in 1:n_plans) { # loop over number conserved
#'   w[[i]] <- list()
#'   for(j in 1:n_trials) { # loop over trials
#'     w[[i]][[j]] <- matrix(rep(625, 16), nrow = 1)
#'     w[[i]][[j]][-sample(1:16, num_pops[i])] <- 5
#'   }
#'  }
#' arma_env_params <- list(mean_value = 16, ar = 0.1, sigma_env = 2, ma = 0)
#'
#' x_arma_sp <- run_cons_plans(w, env_type = "arma", env_params = arma_env_params)
#'
#' plot_cons_plans(x_arma_sp$plans_mv, plans_name = plans_name_sp, cols =
#'   cols, add_all_efs = FALSE, xlim = c(0.02, 0.15), ylim = c(-0.017,
#'     0.017), add_legend = FALSE)
#'
#' # In this version, the pops are wiped out; total abundance changes
#' n_trials <- 50 # number of trials at each n conservation plan
#' num_pops <- c(2, 4, 8, 16) # n pops to conserve
#' n_plans <- length(num_pops) # number of plans
#' w <- list()
#' for(i in 1:n_plans) { # loop over number conserved
#'  w[[i]] <- list()
#'  for(j in 1:n_trials) { # loop over trials
#'    w[[i]][[j]] <- matrix(rep(1000, 16), nrow = 1)
#'    w[[i]][[j]][-sample(1:16, num_pops[i])] <- 5
#'  }
#' }
#' plans_name_n <- paste(num_pops, "populations")
#' arma_env_params <- list(mean_value = 16, ar = 0.1, sigma_env = 2, ma = 0)
#'
#' x_arma_n <- run_cons_plans(w, env_type = "arma", env_params =
#'   arma_env_params, max_a = thermal_integration(16))
#'
#' plot_cons_plans(x_arma_n$plans_mv, plans_name = plans_name_n, cols =
#'   cols, add_all_efs = FALSE, xlim = c(0.02, 0.15), ylim = c(-0.017,
#'     0.017), add_legend = FALSE)
#' }

run_cons_plans <- function(w, env_type, env_params, show_progress =
  TRUE, burn = 1:30, assess_freq = 5, risk_fn = var, ...) {

  plans_mv_n <- list()
  plans_port_n <- list()
  for(i in 1:length(w)) {
    plans_mv_n[[i]] <- list()
    plans_port_n[[i]] <- list()
    n_trials <- length(w[[i]])
    for(j in 1:n_trials) {
      temp <- get_conserv_plans_mv(
                    weights       = w[[i]][[j]],
                    reps          = 1,
                    env_type      = env_type,
                    env_params    = env_params,
                    burn          = burn,
                    assess_freq   = assess_freq,
                    risk_fn       = risk_fn, ...
                            )
    plans_mv_n[[i]][[j]] <- temp$port_mv
    plans_port_n[[i]][[j]] <- temp$port_out
    }
    if(show_progress)
      print(paste("Completed", i, "of", length(w), "conservation plans to evaluate"))
  }

  # here we have a nested list that is 4 elements deep
  # e.g. plans_port_n[[1]][[1]][[1]][[1]]$A
  # let's make that slightly more sane by removing the extra two
  # levels at the end:
  # So, the structure will be plans_port_n[[w]][[trial]]
  plans_port_n_out <- list()
  for(i in 1:length(w)) {
    plans_port_n_out[[i]] <- list()
    for(j in 1:n_trials) {
      plans_port_n_out[[i]][[j]] <- plans_port_n[[i]][[j]][[1]][[1]]
    }
  }

  # Move to list of dataframes instead of list of list of dataframes
  # this is to match the typical output with even population numbers
  # across runs.
  plans_mv_n_dfs <- list()
  for(i in 1:length(w)) {
    plans_mv_n_dfs[[i]] <- list()
    for(j in 1:n_trials) {
      plans_mv_n_dfs[[i]][[j]] <- plans_mv_n[[i]][[j]][[1]]
    }
  }
  for(i in 1:length(w)) {
    plans_mv_n_dfs[[i]] <- do.call("rbind", plans_mv_n_dfs[[i]])
  }

  return(list(plans_mv = plans_mv_n_dfs, plans_port =
      plans_port_n_out))
}
