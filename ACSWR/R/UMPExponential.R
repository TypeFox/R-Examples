UMPExponential <-
function(theta0, n, alpha){
  t <- qgamma(1-alpha, shape=n,scale=theta0)
  return(t)
}
