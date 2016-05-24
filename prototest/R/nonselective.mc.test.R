#### perfroms the non-selective MC test
#### assumes that y is a MATRIX, one column per replication
#### test.index is the index into x for which we test the nullity of the regression coefficient
#### x is assumed to have at most one column from each of the original groups
#### sigma is assumed known -- so a normal reference distribution
nonselective.mc.test <-
function(x, y, test.index, mu, sigma){
  n = nrow (y)
  
  
  ### are we estimating mu?
  X = x
  if (is.null(mu)){
    X = cbind (1, x)
    test.index = test.index+ifelse(is.null(mu), 1, 0) # move the index over one to account for the intercept
  }else{
    y = y - mu
  }
  M = ncol(X)
  
  ### coefficient and standard error
  XX = t(X)%*%X
  Xy = t(X)%*%y
  beta.hat = solve(XX, Xy)[test.index,]
  std.error = sqrt(solve (XX, sigma^2*diag(M))[test.index, test.index])
  
  ### test stat and p.val
  ts = beta.hat/std.error
  p.val = 2*pnorm(abs(ts), 0, 1, lower.tail=FALSE)
  
  list (ts=ts, p.val=p.val)
}
