#### computes the approximate likelihood for given theta, my, sigma
#### data in via y, P and M
compute.approx.lr <-
function(theta, mu, sigma, y, Py, M){
  ### precompute some things
  n = length(y)
  
  -M*theta - 0.5*M*theta^2 - n*log(sigma) - 0.5*sum((y -theta*Py - mu)^2)/sigma^2
}
