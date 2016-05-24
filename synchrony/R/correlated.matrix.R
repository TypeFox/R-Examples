correlated.matrix <- function (rho = 0, sigma = 1, mu = 0, ntimes = 200, nspecies = 10) {
  
  if (length(rho) > 1) {
    corr.mat=matrix(NA, nrow=nspecies, ncol=nspecies)
    corr.mat[upper.tri(corr.mat)]=rho
    corr.mat[lower.tri(corr.mat)]=corr.mat[upper.tri(corr.mat)]
  }
  else {
    corr.mat=matrix(rho, nrow=nspecies, ncol=nspecies)
  }
  diag(corr.mat)=1
  # Cholesky decomposition
  L=try(chol(corr.mat), silent=TRUE)
  if (class(L)!="try-error") {
    community=matrix(rnorm(ntimes*nspecies, sd=1, mean=0), nrow=ntimes, ncol=nspecies) %*% L
    community=scale(community, center=TRUE, scale=TRUE)*sigma+mu
    attr(community, "scaled:center")=NULL
    attr(community, "scaled:scale")=NULL
  }
  else {
    community=NA
    warning("Unable to generate desired correlation matrix (leading minor is not positive definite)")
  }
  
  results=list(rho=rho, sigma=sigma, mu=mu, community=community)
  class(results)="cormat"
  return (results)
}
