`morisita.horn` <-
function(x,y) {
  aN <- sum(x)
  bN <- sum(y)
  da <- sum(x^2)/aN^2
  db <- sum(y^2)/bN^2
  top <- sum(x*y)
  mor <- 2*top/((da+db)*aN*bN)
  return(mor)
}

