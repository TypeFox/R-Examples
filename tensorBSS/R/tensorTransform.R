tensorTransform <-
function(x, A, m){
  r <- length(dim(x)) - 1
  x <- tensor(x, A, alongA = m, alongB = 2)
  aperm(x, cyclicalPermute(m, r))
}
