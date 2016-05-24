update.sigma <-
function(theta, mu, y, Py){
  sqrt(mean((y - theta*Py - mu)^2))
}
