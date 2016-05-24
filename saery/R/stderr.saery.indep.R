stderr.saery.indep <-
function(fit) {
  Finv <- solve(fit[[2]])
  sigma1.std.err <- sqrt(Finv[1,1])
  sigma2.std.err <- sqrt(Finv[2,2])
  beta.std.err <- sqrt(as.vector(diag(fit[[5]])))
  
  return( list(beta.std.err, c(sigma1.std.err, sigma2.std.err)) )
  
}
