#### computes the exact likelihood for given theta, my, sigma
#### data in via y, P and M
compute.exact.lr <-
function(theta, mu, sigma, y, Py, M){
  ### precompute some things
  n = length(y)
  
  M*log(1-theta) - n*log(sigma) - 0.5*sum((y -theta*Py - mu)^2)/sigma^2
}
