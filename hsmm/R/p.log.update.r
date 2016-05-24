#####################################################################################
# Update functions for the parameters which for an explicit solution does not exist #
# Solutions are calculated inside EM using root-finding algorithms                  #
#####################################################################################

# Update function for p.log
# ------------------------
p.log.update <- function(p, j, M, eta){
  x <- sum(eta[1:M, j] * (1 / p * 1:M + 1 / (1 - p) / log(1 - p)))
  return(x)
  }
