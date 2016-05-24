parse.joint <-
function(P, v, cond, var) {
  P.num <- P
  P.num$sumset <- c(union(P$sumset, setdiff(var, union(v, cond))))
  if (length(cond) > 0) {
    P.den <- P
    P.den$sumset <- c(union(P$sumset, setdiff(var, cond)))
    if (length(P.den$children) > 0) {
      P.num$fraction <- TRUE
      P.num$divisor <- P.den
    }
  }
  return(P.num)
}
