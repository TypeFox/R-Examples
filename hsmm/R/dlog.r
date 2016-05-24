#######################################################
# Density functions of the sojourn time distributions #
#######################################################

# Density function of the logarithmic distribution
# ------------------------------------------------
dlog <- function(x, p){
  y <- - p^(x + 1) / ((x + 1) * log(1 - p))
  return(y)
  }
