cumDistFun = function(X, x, px){
  Px = cumsum(px)
  r = rep(0, length(X))
  r[X >= max(x)] = 1
  i = X >= min(x) & X < max(x)
  j = sapply(X[i], function(y)max(which(x <= y)))
  r[i] = Px[j]
  return(r)
}