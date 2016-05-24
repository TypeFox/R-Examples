### Generate phi.Obs based on param returned from mixnorm.optim().

### log(Phi) ~ \sum_k p_k N(mu_k, sigma_k^2)
### log(Phi^{obs}) = log(Phi) + N(0, sigma_e^2)
simu.mixnormerr <- function(n, param){
  id.K <- sample(1:param$K, n, replace = TRUE, prob = param$prop)
  n.K <- tabulate(id.K, nbins = param$K)

  ret <- list(Phi = rep(0.0, n), phi.Obs = rep(0.0, n), id.K = id.K)
  for(i.k in 1:length(n.K)){
    if(n.K[i.k] != 0){
      tmp <- rnorm(n.K[i.k], param$mu[i.k], sqrt(param$sigma2[i.k]))
      tmp.err <- tmp + rnorm(n.K[i.k], 0, sqrt(param$sigma2.e))

      ret$Phi[id.K == i.k] <- exp(tmp)
      ret$phi.Obs[id.K == i.k] <- exp(tmp.err)
    }
  }

  ret
} # End of simu.mixnormerr().
