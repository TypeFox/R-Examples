LikFunk <- function(W, psi, dX){

# Computes value of log-likelihood at psi where
#
# L(psi) := sum(W * psi) - sum_{X_1} ^ {X_n} exp(psi(z))
#         = sum(W * psi) - sum_{k = 1} ^ {m - 1} J(psi_k, psi_(k+1))
#           + sum_{k = 2} ^ {m - 1} exp(psi_k)

m <- length(psi)
last <- 0
if (m > 2){last <- sum(exp(psi[2:(m - 1)]))}
L <- sum(W * psi) - sum(J00(psi[1:(m - 1)], psi[2:m], dX)) + last

return(L)
}
