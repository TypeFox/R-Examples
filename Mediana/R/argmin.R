######################################################################################################################

# Function: argmin.
# Argument: p, Vector of p-values (1 x m)
#           w, Vector of hypothesis weights (1 x m)
#           processed, Vector of binary indicators (1 x m) [1 if processed and 0 otherwise].
# Description: Hidden function used in the Chain function. Find the index of the smallest weighted p-value among the non-processed null hypotheses with a positive weight (index=0 if
# the smallest weighted p-value does not exist) in a chain procedure

argmin = function(p, w, processed) {

  index = 0
  m = length(p)
  for (i in 1:m) {
    if (w[i] > 0 & processed[i] == 0) {
      if (index == 0) {
        pmin = p[i]/w[i]
        index = i
      }
      if (index > 0 & p[i]/w[i] < pmin) {
        pmin = p[i]/w[i]
        index = i
      }
    }
  }
  return(index)
}
# End of argmin