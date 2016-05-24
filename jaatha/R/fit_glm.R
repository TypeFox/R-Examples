fit_glm <- function(x, sim_data, ...) UseMethod("fit_glm")
fit_glm.default <- function(x, sim_data, ...) {
  stop("Unknown Summary Statistic")
}


#' @export
fit_glm.jaatha_model <- function(x, sim_data, ...) { #nolint
  "Fits a GLM to the simulation results"
  lapply(x$get_sum_stats(), fit_glm, sim_data, ...)
}


#' @export
fit_glm.jaatha_stat_basic <- function(x, sim_data, ...) {
  "Fits a GLM for each entry of the simulation results"
  Y <- do.call(rbind, lapply(sim_data, function(data) data[[x$get_name()]]))
  X <- cbind(1, 
             do.call(rbind, lapply(sim_data, function(data) data$pars_normal)))
  
  glms <- lapply(1:ncol(Y), function(i) {
    suppressWarnings(
      stats::glm.fit(X, Y[, i], family = stats::poisson("log"), 
                     control = list(maxit = 100))[c("coefficients", 
                                                    "converged")]
    )
  })
  
  vapply(glms, function(x) {
    if (!x$converged) stop("GLM did not converge", call. = FALSE)
    numeric(0)
  }, numeric(0))
  
  glms
}
