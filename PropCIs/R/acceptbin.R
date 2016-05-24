acceptbin <-
function(x, n, p){

  p1 = 1 - pbinom(x - 1, n, p)
  p2 = pbinom(x, n, p)
  a1 = p1 + pbinom(qbinom(p1, n, p) - 1, n, p)
  a2 = p2 + 1 - pbinom(qbinom(1 - p2, n, p), n, p)
  return(min(a1,a2))
}

