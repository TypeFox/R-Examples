# Special linear time order function, works only when x
# is a scrambled vector of integers.

Order <- function(x) {
  n = length(x)
  o = numeric(n)
  o[x] = Seq(1,n)
  return(o)
}
