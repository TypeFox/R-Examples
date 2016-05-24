################################################################################
#
#  Generate saved objects for slp package
#
#  This routine is never intended to be run by a user. It generates the saved
#  basis objects by running slp() for certain specific N, W and K. 
#
#  ** Note: W is stored as df/year, integer
#
#  (C) wburr, July 2014
#  Licensed under GPL-2
#
################################################################################

.generateSaved <- function(fName, workDir) {

  setwd(workDir)

  # This lists all the saved objects by N, W and two choices of K: 2NW - 1 and 2NW (approx, rounded to integers)
  N = c(365, 366, 730, 731) 
  W = c(6, 7)

  nSaved <- 2 * length(N) * length(W)

  slpSavedObjects <- vector("list", nSaved)

  idx <- 1
  for(j in 1:length(N)) {
      for(k in 1:length(W)) {
          slpSavedObjects[[idx]]   <- c(N[j], W[k], floor(2 * N[j] * W[k] / 365.2425))
          slpSavedObjects[[idx+1]] <- c(N[j], W[k], ceiling(2 * N[j] * W[k] / 365.2425))
          idx <- idx + 2
      }
  }
  save(file = fName, slpSavedObjects)

  # Now generate these objects (WARNING: THIS TAKES A _VERY_ LONG TIME)

  for(j in 1:length(N)) {
      for(k in 1:length(W)) {
         cat(paste(j, " - ", k, "\n"))

         K <- floor(2 * N[j] * W[k] / 365.2425)
         basis <- slp(1:N[j], W = W[k] / 365.2425, K = K, intercept = TRUE, forceC = TRUE) 
         save(file = paste0("basis_N_", N[j], "_W_", W[k], "_K_", K, ".RData"), basis)

         K <- ceiling(2 * N[j] * W[k] / 365.2425)
         basis <- slp(1:N[j], W = W[k] / 365.2425, K = K, intercept = TRUE, forceC = TRUE) 
         save(file = paste0("basis_N_", N[j], "_W_", W[k], "_K_", K, ".RData"), basis)
      }
  }
}
