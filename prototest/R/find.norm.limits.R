#### computes the limits of the norm of Pz = P(y - mu)
#### used for non-sampling selective inference when mu is specified
#### selective constarints such that Ay <= b
#### y is a VECTOR
find.norm.limits <-
function(U, y, mu, A, b){
  # precompute
  b.tilde = b - mu*apply(A, 1, sum)
  z = y - mu
  Pz = U%*%(t(U)%*%z)
  APz = A%*%Pz
  delta = A%*%z - APz
  v = Pz/sqrt(sum(Pz^2))
  
  Av = A%*%v
  
  ## limits for each of the replications of y
  # positive and negtaive indices
  which.pos = which (Av > 0)
  which.neg = which (Av < 0)
  values = (b.tilde - delta)/Av
  
  # compute the limits
  options(warn=-1)
  lower = max(values[which.neg])
  upper = min(values[which.pos])
  options(warn=0)
  
  c(lower, upper)
}
