Seq <- function(a, b, ...) {
  if (a<=b) return(seq(a,b,...))
  else return(numeric(0))
}
