#' Bayesian inference on a mutlivariate normal (MVN) mean with a multivariate normal (MVN) prior
#' 
#' Evaluates posterior density for \eqn{\mu}{mu}, the mean of a
#' MVN distribution, with a MVN prior on \eqn{\mu}{mu}
#' 
#' 
#' @param y a vector of observations from a MVN distribution with unknown
#' mean and known variance-covariance.
#' @param m0 the mean vector of the MVN prior, or a scalar constant so that the prior
#' vector of length \eqn{k}{k} with the same element repeated k times, e.g. \code{m0 = 0}
#' @param V0 the variance-covariance matrix of the MVN prior, or the diagonal 
#' of the variance-covariance matrix of the MVN prior, or a scalar constant, say \eqn{n_0}{n0},  
#' so the prior is \eqn{n_0\times \mathbf{I}_k}{n0 * I} where \eqn{\mathbf{I}_k}{I} is the \eqn{k}{k} by \eqn{k}{k} identity matrix.
#' @param Sigma the known variance covariance matrix of the data. If
#' this value is NULL, which it is by default, then the sample covariance is used. NOTE:
#'  if this is the case then the cdf and quantile functions should really be multivariate
#'  t, but they are not - in which case the results are only (approximately) valid for large samples. 
#' @return A list will be returned with the following components: 
#' \item{mean}{the posterior mean of the MVN posterior distribution}
#' \item{var}{the posterior variance-covariance matrix of the MVN posterior distribution}
#' \item{cdf}{a function that will evaluation the posterior cdf at a given point. This function calls \code{mvtnmorm::pmvnorm}.}
#' \item{quantile}{a function that will find quantiles from the posterior given input probabilities. This function calls \code{mvtnorm::qmvnorm}.}
#' @keywords misc
#' @export 
mvnmvnp = function(y, m0 = 0, V0 = 1, Sigma = NULL){
  
  yBar = matrix(colMeans(y), ncol = 1)
  k = ncol(y)
  
  if(length(m0) == 1){
    m0 = rep(m0, k)
  }
  m0 = matrix(m0, ncol = 1)
  
  if(!is.matrix(V0)){
    if(length(V0) == 1){
      V0 = diag(V0, k, k)
    }else{
      if(length(V0) != k){
        stop("V0 must either be a scalar, a vector of length k, or a k x k symmetric matrix, where k = ncol(y)")
      }else{
        V0 = diag(V0, k, k)
      }
    }
  }else{
    if(!isSymmetric(V0)){
      stop("The prior variance-covariance V0 must be a k x k symmetric matrix, where k = ncol(y)")
    }
  }
  
  cat("The prior mean is:\n\n")
  cat(paste0(m0, collapse = " "))
  cat("\n\n")
  cat("The prior variance is:\n\n")
  print(V0)
  cat("\n\n")
  
  
  n = nrow(y)
  
  if(is.null(Sigma)){
    Sigma = cov(y)
    cat("Using the sample variance for Sigma\n")
  }
  
  Sigma.inv = solve(Sigma)
  
  prior.precision = solve(V0)
  post.precision = prior.precision + n * Sigma.inv
  post.var = solve(post.precision)
  post.mean = post.var %*% prior.precision %*% m0 + n * post.var %*% Sigma.inv %*% yBar
  
  cat("The posterior mean is:\n\n")
  cat(paste0(post.mean, collapse = " "))
  cat("\n\n")
  cat("The posterior variance is:\n\n")
  print(post.var)
  
  results = list(parameter = 'mu',
                 mean = post.mean, 
                 var = post.var,
                 cdf = function(x, ...)mvtnorm::pmvnorm(x, post.mean, post.var, ...),
                 quantileFun = function(probs, ...)mvtnorm::qmvnorm(probs, post.mean, post.var, ...))
  
  class(results) = 'Bolstad'
  
  invisible(results)
}
