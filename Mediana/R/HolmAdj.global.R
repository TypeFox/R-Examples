######################################################################################################################

# Function: HolmAdj.global.
# Argument: p, Vector of p-values (1 x m)
#           n, Total number of testable hypotheses (in the case of modified mixture procedure) (1 x 1)
#       gamma, Vector of truncation parameter (1 x 1)
# Description: Compute global p-value for the truncated Holm multiple testing procedure. The function returns the global adjusted pvalue (1 x 1)

HolmAdj.global = function(p, n, gamma) {

  # Number of p-values
  k = length(p)
  if (k > 0 & n > 0) {
    if (gamma == 0)
    {
      adjp = n * min(p)
    }  # Bonferonni procedure
    else if (gamma <= 1) {
      # Truncated Holm procedure Index of ordered pvalue
      ind = order(p)
      # Denominator (1 x m)
      seq = seq_vector(k)
      denom = gamma/(k - seq + 1) + (1 - gamma)/n
      # Adjusted p-values
      sortp = sort(p)
      adjp = min(cummax(sortp/denom)[order(ind)])
    }
  } else adjp = 1

  return(adjp)

}
# End of HolmAdj.global