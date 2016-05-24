update.theta <-
function(mu, sigma, yPy, yP1, M){
  1 - mu*yP1/2/yPy - sqrt(mu^2*yP1^2 + 4*sigma^2*M*yPy)/2/yPy
}
