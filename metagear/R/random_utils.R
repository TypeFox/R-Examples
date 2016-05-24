
# Function: GenerateMultivariatePoisson
# Yahav, I. and G. Shmueli (2011) On generating multivariate Poisson data
# in management science applications. Appl. Stochastic Models Bus. Ind. 28: 91-102. 
GenerateMultivariatePoisson <- function(p, samples, R, lambda){
  normal_mu=rep(0, p)
  normal = mvrnorm(samples, normal_mu, R)
  unif=pnorm(normal)
  pois=t(qpois(t(unif), lambda))
  return(pois)
}
