graphlaplacian <-
function(A, zeta=0.01){
  if (sum(A != t(A)) > 0){
    stop("This method only works for symmetric A!")
  }
  Adeg = apply(abs(A), 1, sum)   # a p-dim vector
  Adeg[Adeg==0] = 1
  AdegInv = (Adeg + zeta)^(-0.5)
  Lunnorm = diag(Adeg) - A
  Lnorm = diag(AdegInv) %*% A %*% diag(AdegInv)
  
  return(list(Lunnorm = Lunnorm, Lnorm = Lnorm))
}
