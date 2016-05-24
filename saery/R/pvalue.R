pvalue <-
function(beta.fitted, fit) {
  z <- abs(beta.fitted)/sqrt(as.vector(diag(fit[[5]])))
  p <- pnorm(z, lower.tail=F)
  
  return( 2*p )
  
}
