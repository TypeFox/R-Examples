mc.selection.A.b <-
function(y, X, mu){
  # extract some parameters
  n = length(y)
  p = ncol(X)
  
  # standardise X
  X.tilde = apply (X, 2, function(col){col-mean(col)})
  X.tilde = apply (X.tilde, 2, function(col){col/sqrt(sum(col^2))})
  
  # is mu specified? adjust y accordingly
  if (is.null(mu)){
    y.tilde = y - mean(y)
  }else{
    y.tilde = y - mu
  }
  
  ### find the maximum absolute correlation
  which.col = which.max (abs(t(X.tilde)%*%y.tilde))
  max.val = sum(X.tilde[,which.col]*y.tilde)
  max.sign = sign (max.val)
  
  ### A, b matrices
  x.i = X.tilde[,which.col,drop=FALSE]
  X.min.i = X.tilde[,-which.col, drop=FALSE]
  ones = rep (1, p-1)
  
  A = rbind (
    t(X.min.i) - max.sign*(ones%*%t(x.i)),
    -t(X.min.i) - max.sign*(ones%*%t(x.i)),
    -max.sign*t(x.i)
  )
  b = rep (0, nrow(A))
  
  ## adjust A and b according to the specification of mu
  if (is.null(mu)){ # A changes; not b
    A = A - apply(A, 1, sum)%*%t(rep(1, ncol(A)))/n
  }else{
    b = b + mu*apply(A, 1, sum)
  }
  
  list(which.col=which.col, max.val=max.val, A=A, b=b)
}
