comb.summands <-
function(n, k) {
# Returns all combinations of k numbers between 0 and n whose sum equals n.
#
# This function was written by Christoph Kurz.
  recursion <- function(n, k) {
    if (k == 1)
      list(n)
    else
      unlist(lapply(0:n, function(i) Map(c, i, recursion(n - i, k - 1))), recursive = F)
  }
    matrix(unlist(recursion(n, k)), byrow = T, ncol = k)
}
