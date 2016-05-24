######################################################################################################################

# Function: BonferroniAdj.global.
# Argument: p, Vector of p-values (1 x m)
#           n, Total number of testable hypotheses (in the case of modified mixture procedure) (1 x 1)
#       gamma, Vector of truncation parameter (1 x 1)
# Description: Compute global p-value for the Bonferroni multiple testing procedure. The function returns the global adjusted pvalue (1 x 1)

BonferroniAdj.global = function(p, n, gamma) {

  # Number of p-values
  k = length(p)
  if (k > 0 & n > 0) {
    adjp = n * min(p)  # Bonferonni procedure
  } else adjp = 1
  return(adjp)
}
# End of BonferroniAdj.global