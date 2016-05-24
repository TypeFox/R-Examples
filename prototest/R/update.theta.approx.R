update.theta.approx <-
function(mu, sigma, yPy, yP1, M){
  (yPy - yP1*mu - sigma^2*M)/(yPy + sigma^2*M)
}
