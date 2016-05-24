######################################################################################################################

# Function: HommelAdj.global.
# Argument: p, Vector of p-values (1 x m)
#           n, Total number of testable hypotheses (in the case of modified mixture procedure) (1 x 1)
#       gamma, Vector of truncation parameter (1 x 1)
# Description: Compute global p-value for the truncated Hommel multiple testing procedure. The function returns the global adjusted pvalue (1 x 1)

HommelAdj.global = function(p, n, gamma) {

  # Number of p-values
  k = length(p)
  if (k > 0 & n > 0) {
    if (gamma == 0)
    {
      adjp = n * min(p)
    }  # Bonferonni procedure
    else if (gamma <= 1) {
      # Truncated Hommel procedure
      seq = 1:k
      denom = seq * gamma/k + (1 - gamma)/n
      sortp = sort(p)
      adjp = min(sortp/denom)
    }
  } else adjp = 1
  return(adjp)
}
# End of HommelAdj.global