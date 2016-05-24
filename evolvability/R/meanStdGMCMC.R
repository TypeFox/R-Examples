#' @export
meanStdGMCMC = function(G_mcmc, means_mcmc){
  X = cbind(means_mcmc, G_mcmc)
  n1 = ncol(means_mcmc)
  n2 = ncol(X)
  t(apply(X, 1, function(x) x[(n1+1):n2]/c(x[1:n1]%*%x[1:n1])))
}
