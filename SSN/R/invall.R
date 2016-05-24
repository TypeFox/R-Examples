invall <-
function(X) {
  indmginv <- min(svd(X)$d) < 1e-10
  if(indmginv == TRUE) mginv(X)
  else(solve(X))
}

