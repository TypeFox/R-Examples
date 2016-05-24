#####################################################################################
# Update functions for the parameters which for an explicit solution does not exist #
# Solutions are calculated inside EM using root-finding algorithms                  #
#####################################################################################

# Update function for nu
# ------------------------
nu.update <- function(x, dblPara){
  x <- -digamma(exp(x) / 2) + log(exp(x) / 2) + 1 + dblPara + digamma((exp(x) + 1) / 2) - log((exp(x) + 1) / 2);
  return(x)
  }
