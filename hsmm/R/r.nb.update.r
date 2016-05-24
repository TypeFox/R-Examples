#####################################################################################
# Update functions for the parameters which for an explicit solution does not exist #
# Solutions are calculated inside EM using root-finding algorithms                  #
#####################################################################################

# Update function for r.nb
# ------------------------
r.nb.update <- function(r, j, M, eta){
  y <- sum(eta[1:M, j])
  z <- sum(eta[1:M, j] * (1:M - 1 + r))
  x <- sum(eta[1:M, j] * (digamma(1:M + r - 1) - digamma(r) + log(r * y / z)))
  return(x)
  }
