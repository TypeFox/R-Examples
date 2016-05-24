#' Get coda object
#'
#' Returns an \code{mcmc} object for use in \code{coda} functions for 
#' convergence diagnostics, HPD intervals, etc.
#'
#' @param model a \code{bfa} model
#' @param loadings Return samples of the factor loadings? (default is TRUE)
#' @param scores Return samples of the factor scores? (default is FALSE)
#' @param scale Return factor loadings on the correlation scale? (default is TRUE)
#' @param positive Post-process to enforce positivity constraints
#' @return An \code{mcmc} object
#' @export

get_coda <- function(model, loadings=TRUE, scores=FALSE, scale=TRUE, positive=FALSE) {
  k = model$K; p = model$P; n=model$N
  nsim = model$nsim; thin = model$thin; nburn = model$nburn
  samples = NULL
  kn=c()
  for (i in 1:k) kn = c(kn, rep(i, p))
  loadnames = paste(model$varlabel, kn)
  
  if (loadings) {
    pl = model$post.loadings
    if (scale) {
      if (dim(pl)[2]>1) {
        for (i in 1:dim(pl)[3]) {
          pl[,,i] = pl[,,i]/sqrt(model$post.sigma2[i,]+rowSums(pl[,,i]^2))
        }
      } else {
        for (i in 1:dim(pl)[3]) {
          pl[,,i] = pl[,,i]/sqrt(model$post.sigma2[i,]+pl[,,i]^2)
        }
      }
      if (positive) {
        for (i in 1:k) {
          row = which(model$loadings.restrict[,k]==2)
          ind = which(model$post.loadings[row,k,]<0)
          model$post.loadings[,k,ind] = -1*model$post.loadings[,k,ind]
        }
      }
    }
    dim(pl) = c(1,p*k, nsim/thin)
    pl = pl[1,,]
    rownames(pl) = loadnames
    samples = cbind(samples, t(pl))
  }
  
  kn=c()
  for (i in 1:k) kn = c(kn, rep(i, n))
  scorenames = paste(model$obslabel, kn)
  if(scores) {
    ps = model$post.scores
    dim(ps) = c(1,n*k, nsim/thin)
    ps = ps[1,,]
    rownames(ps) = scorenames
    samples = cbind(samples, t(ps))
  }
  
  codaobj = mcmc(samples, nburn+1, nsim+nburn, thin=thin)
}