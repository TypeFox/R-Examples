# This is a hidden function of the l2boost package.

# l2boost function to determine step size for long descent from Theorem 2 of Ehrlinger and Ishwaran (2012)
# |Dk| must be less than or equal to 1; We force it if otherwise, to enforce step=infty when dk=Rk=1
#
# @param rho.m vector of stagewise regression parameters 
# @param Rk correlation coefficient of candidate and current directions
# @param lr current direction index
# @param nu l1 shrinkage parameter
# 
# 
# 
# @references John Ehrlinger, Hemant Ishwaran (2012). Characterizing l2boosting. \emph{Annals of Statistics}, 40 (2), 1074-1101

mstep.long <- function(rho.m, Rk, lr, nu) {
    Dk <- pmin(pmax(rho.m/rho.m[lr], -1, na.rm = TRUE), 1)
    num <- Dk - Rk
    nu.r <- abs(num)/(1 - Rk * sign(num))
    mjk <- floor(1 + log(nu.r)/log(1 - nu))
    mjk[Dk == 1 & Rk == 1] <- Inf
    mjk
}
