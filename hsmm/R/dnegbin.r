#######################################################
# Density functions of the sojourn time distributions #
#######################################################

# Density function of the neg. bin. distribution
# ----------------------------------------------
dnegbin <- function(x, r, varpi){
  y <- exp(lgamma(x + r) - lgamma(r) - lgamma(x + 1)) * varpi^r * (1 - varpi)^x
  return(y)
  }
