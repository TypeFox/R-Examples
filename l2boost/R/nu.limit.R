# This is a hidden function of the l2boost package.

# limiting nu size, this comes directly from the RHS of theorem 2 of Ehrlinger and Ishwaran 2012.
#
# @param rho.m vector of stagewise regression parameters
# @param Rk correlation coefficient of current and candidate coordinate directions
# @param lr candidate direction index
# 
# @references Ehrlinger and Ishwaran (2012). Characterizing l2boosting. \emph{Annals of Statistics}, 40 (2), 1074-1101
nu.limit <- function(rho.m, Rk, lr) {
  Dk <- rho.m/rho.m[lr]
  num <- Dk - Rk
  return(1 - abs(num)/(1 - Rk * sign(num)))
}
