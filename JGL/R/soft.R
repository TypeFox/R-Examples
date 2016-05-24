soft <-
function(a,lam,penalize.diagonal){ # if last argument is FALSE, soft-threshold a matrix but don't penalize the diagonal
  out <- sign(a)*pmax(0, abs(a)-lam)
  if(!penalize.diagonal) diag(out) <- diag(a)
  return(out)
}

