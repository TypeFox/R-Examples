stderr.saery.MA1 <-
function(fit) {
  Finv <- solve(fit[[2]])
  sigma1.std.err <- sqrt(Finv[1,1])
  sigma2.std.err <- sqrt(Finv[2,2])
  theta.std.err <- sqrt(Finv[3,3])
  beta.std.err <- sqrt(as.vector(diag(fit[[5]])))
  
  return( list(beta.std.err, c(sigma1.std.err, sigma2.std.err, theta.std.err)) )
  
}
