# Make sure that no two rows of A are the same (this
# works with probability one).

checkrows <- function(A) {
  b = rnorm(ncol(A))
  a = sort(A%*%b)
  if (any(diff(a)==0)) return(TRUE)
  return(FALSE)
}
