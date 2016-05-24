potence <-
function(a, b) {
  c = a
  if(b < 0 || a < 0)
    stop("invalid input")
    if( b == 0) return(1)
      while(b > 1) {
        c = c*a
        b = b-1
      }
  return(c)
}
