### finds the limits of the approx LR stat (not the statistic itself, but the chi-squared part of it)
### given limits for the norm of Pz, as found by find.norm.limits
limits.approx.lr <-
function(U, y, A, b, mu, sigma){
  ### norm limits
  norm.limits = find.norm.limits (U=U, y=y, mu=mu, A=A, b=b)
  
  norm.limits^2/sigma^2
}
