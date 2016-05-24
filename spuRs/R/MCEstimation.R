# program spuRs/resources/scripts/MCEstimation.R
# loadable spuRs function

# This function estimates the transition matrix for a discrete MC
# with known state space {0,1,..,n}
# based on observations of the chain stored in statehist.

MCEstimation <- function(statehist, n) {
  ntrans <- length(statehist) - 1
  TC <- matrix(0, nrow = n+1, ncol = n+1)
  EP <- TC
  for (trans in 2:ntrans) {
    TC[statehist[trans - 1] + 1, statehist[trans] + 1] <-
      TC[statehist[trans - 1] + 1, statehist[trans] + 1] + 1
  }
  for (row in 1:(n+1)) {
    EP[row,] <- TC[row,]/rowSums(TC)[row]
  }
  return(EP)
}