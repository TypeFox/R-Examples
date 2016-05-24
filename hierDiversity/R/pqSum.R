pqSum <-
function(p, q = q) {
  if (q != 1) {
    temp.sum <- sum(p[p > 0]^q)
  }
  if (q == 1) {
    temp.sum <- sum(p[p > 0] * log(p[p > 0]))
  }
  return(temp.sum)
}
