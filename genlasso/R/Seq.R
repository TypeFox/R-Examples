# Returns a sequence of integers from a to b if a <= b,
# otherwise nothing. You have no idea how important this
# function is...

Seq <- function(a, b, ...) {
  if (a<=b) return(seq(a,b,...))
  else return(numeric(0))
}
