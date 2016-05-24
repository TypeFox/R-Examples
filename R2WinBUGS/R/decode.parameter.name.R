"decode.parameter.name" <-
function (a){
#
# Decodes Bugs parameter names
#   (e.g., "beta[3,14]" becomes "beta" with 2 indexes:  3 and 14)
# for use by the bugs.sim() function
#
  left.bracket <- regexpr ("[[]", a)
  if (left.bracket==-1){
    root <- a
    dimension <- 0
    indexes <- NA
  }
  else {
    root <- substring (a, 1, left.bracket-1)
    right.bracket <- regexpr ("[]]", a)
    a <- substring (a, left.bracket+1, right.bracket-1)
    indexes <- as.numeric(unlist(strsplit(a, ",")))
    dimension <- length(indexes)
  }
  list(root=root, dimension=dimension, indexes=indexes)
}
