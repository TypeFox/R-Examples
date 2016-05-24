### finds the limits of the exact LR stat
### given limits for the norm of Pz, as found by find.norm.limits
limits.exact.lr <-
function(U, y, A, b, mu, sigma, M){
  ### norm limits
  norm.limits = find.norm.limits (U=U, y=y, mu=mu, A=A, b=b)
  
  ### boundary values
  boundary.vals = M*log(sigma^2*M) - M + norm.limits^2/sigma^2 - 2*M*log(norm.limits)
  
  ### is sigma*sqrt(M) in the norm limits?
  if (sigma*sqrt(M) <= norm.limits[2] & sigma*sqrt(M) >= norm.limits[1]){
    lower = 0
  }else{
    lower = min(boundary.vals)
  }
  upper = max(boundary.vals)
  c(lower, upper)
}
