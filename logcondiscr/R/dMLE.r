dMLE <- function(X, W, psi_o = NA, prec = 1e-7){

# dMLE computes vector psi s.t. the function L, defined as
#
# L(psi) := sum(W * psi) - sum_{k = 1} ^ {m - 1} J(psi_k, psi_(k + 1)) 
#                     + sum_{k = 2} ^ {m - 1} exp(psi_k)
#
# is maximized.
# 
# INPUT
# - X     : vector (X_1 < ... < X_m) of observations
# - W     : vector of weights that sum to one
# - psi_o : initialization vector (optional)
# 
# OUTPUT
# - psi   : vector that maximizes L below
# - L     : value of L at psi
#
#
# Version 17.10.07, Kathrin Weyermann
# Ported to R by Kaspar Rufibach, October 2010
########################################################################

# initialization
m <- length(X)
dX <- diff(X)
if (identical(psi_o, NA)){psi_o <- log(1 / m[1]) * rep(1, m)}
psi <- psi_o 
L <- LikFunk(W, psi, dX) 
iterDMLE <- 0
stepDirection <- Direction(W, psi, dX)
norm2 <- sum(stepDirection ^ 2)

while ((norm2 > prec) & (iterDMLE < 50)){
    iterDMLE <- iterDMLE + 1
    stepSize <- StepSize(W, psi, dX, stepDirection)
    psi <- psi + stepSize * stepDirection
    L <- LikFunk(W, psi, dX) 
    stepDirection <- Direction(W, psi, dX)
    norm2 <- sum(stepDirection ^ 2)
} # end while

res <- list("psi" = psi, "L" = L)
return(res)
}

