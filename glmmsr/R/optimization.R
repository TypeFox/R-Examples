#' Maximize the approximated log-likelihood
#'
#' Find the approximate maximum likelihood estimate, and estimate the variance
#' of the estimate
#' @param lfun the approximated loglikelihood function
#' @param p_beta the number of covariates
#' @param p_theta the number of random effects
#' @inheritParams glmm
#' @keywords internal
#' @return A list, containing the parameter estimate and variance matrix
optimize_glmm <- function(lfun, p_beta, p_theta, prev_fit = NULL,
                          verbose = 1L){
  p <- p_theta + p_beta
  if(length(prev_fit) > 0){
    mu <- prev_fit$estim
  } else{
    mu <- c(rep(0.5, p_theta), rep(0, p_beta))
  }

  devfun_ext <- function(param, stop_on_error = FALSE){
    result <- tryCatch(-2 * lfun(param),
                       error = function(e) {
                         if(stop_on_error)
                           stop(e)
                         else
                           return(Inf)
                         })

    if(verbose == 1L) {
      count <<- count + 1
       if(count > print_gap) {
         cat(".")
         count <<- 0
       }
    }
    if(verbose == 2L & started) {
      cat(sprintf("%7.4f", param), " : ", sprintf("%7.4f", -result / 2), "\n")
    }
    result
  }
  devfun_std <- function(param_std, stop_on_error = FALSE){
    param <- as.numeric(param_std + mu)
    return(devfun_ext(param, stop_on_error))
  }

  time_threshold <- 0.1
  time_interval <- 1
  print_gap <- 100
  count <- 0
  started <- FALSE

  par0 <- rep(0, p)
  tryCatch(d0 <- devfun_std(par0, stop_on_error = TRUE),
           error = function(e) {
             stop("Could not approximate the likelihood at the starting parameters for optimization, due to ", e,
                  call. = FALSE)
           })
  t0 <- system.time(devfun_std(par0))[[1]]
  print_gap <- floor(time_interval / t0)

  if(verbose > 0L && t0 > time_threshold)
    cat("Approximating the likelihood at each point takes", t0, "seconds. \n")

  started <- TRUE
  if(verbose > 0L)
    cat("Fitting the model.")
  if(verbose == 2L) {
    cat("\n")
  }
  out_std <- optim(par0, devfun_std, hessian=TRUE, method="BFGS",
                   control = list(maxit = 500))
  if(out_std$convergence != 0)
    warning("optim did not converge")
  if(verbose > 0L){
    cat(" done.\n")
  }
  hess <- out_std$hessian
  estim_std <- out_std$par
  estim = as.numeric(estim_std + mu)

  Sigma <- 2 * solve(hess)
  if(min(abs(eigen(Sigma, only.values=TRUE)$values)) < 1e-6) {
    warning("May have problem in convergence: Sigma has a very small eigenvalue")
  }
  result <- list(estim = estim, Sigma = Sigma)
  return(result)
}
