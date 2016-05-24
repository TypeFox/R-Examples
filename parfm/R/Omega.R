################################################################################
#  Omega matrix for Positive Stable                                            #
################################################################################
#                                                                              #
#                                                                              #
#   Date: February 02, 2012                                                    #
#   Last modification on: February 03, 2012                                    #
################################################################################

Omega <- function(D, correct, nu) {
  Omega <- matrix(NA, nrow=D, ncol=D, dimnames=list(q=1:D, m=0:(D-1)))
  
  Omega[, "0"] <- 10^-correct
  
  if(D > 1) {
    diag(Omega)[-1] <- sapply(2:D, function(q) {
      (1 - nu)^(1 - q) *
        prod(10^(-correct / (q - 1)) *
        (q - 1 + nu - seq(from=1, to=q-1, by=1)))
    })
    if(D > 2)
      for (q in 3:D)  #mPrime = m + 1
        Omega[q, 2:(q - 1)] <- sapply(2:(q - 1), function(mPrime) { 
          Omega[q - 1, mPrime] +
            Omega[q - 1, mPrime - 1] * 
            ((q - 1) / (1 - nu) - (q - (mPrime - 1)))
        })
  }
  
  return(Omega)
}