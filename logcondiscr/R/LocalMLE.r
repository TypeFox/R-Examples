LocalMLE <- function(X, W, isKnot, psi_o, prec){

########################################################################
#
# Compute concave/convex knots and auxiliary vector D
# where D(a) := DL(psi, Delta_a), so the directional derivative in direction
# "new knot at X_a".
#
# Delta_a is defined as
#
# Delta_a(X_j) = 0                                      if X_j is outside X_KK[i], X_KK[i + 1] 
# Delta_a(X_j) ~ (X_j - X_KK[i]) / (X_k - X_KK[i])      if X_j in [X_KK[i], X_k]
# Delta_a(X_j) ~ (X_KK[i] - X_j) / (X_KK[i + 1] - X_k)  if X_j in [X_k, X_KK[i + 1]]
#     
# So
#   D(a) = sum_{j = KK[i]} ^ {KK[i + 1]} wtmp(j) * delta_a(X_j)
# where
#   wtmp = W(ind(2:end - 1))
#          - J10( psi(ind(2:end - 1)), psi(ind(1:end - 2)), dX(ind(1:end - 2)))
#          - J10( psi(ind(2:end - 1)), psi(ind(3:end)), dX(ind(2:end - 1)))
#          + exp( psi(ind(2:end - 1)))
#
# Matlab version 20.11.07, Kathrin Weyermann
# Ported to R by Kaspar Rufibach, October 2010
########################################################################

# initialization
res1 <- LocalCoarsen(X, W, isKnot)
X2 <- res1$X2
W2 <- res1$W2
res2 <- dMLE(X2, W2, psi_o[which(isKnot == 1)])
psi2 <- res2$psi
L <- res2$L
psi <- LocalExtend(X, isKnot, X2, psi2)

# determine concave/convex knots
dX <- diff(X)
conc <- isKnot * LocalConcavity(psi, dX)

# auxiliary vector D
D <- rep(0, length(X))
KK <- which(isKnot == 1)

for (i in 1:(length(KK) - 1)){
    if (KK[i + 1] > (KK[i] + 1)){
        dtmp <- X[KK[i + 1]] - X[KK[i]]
        
        # indices between two knots, incl. the knots
        ind <- KK[i]:KK[i + 1]
        xtmp <- (X[ind[2:(length(ind) - 1)]] - X[KK[i]]) / dtmp
        tmp_grad <- GradientL(W[ind], psi[ind], dX[ind[1:(length(ind) - 1)]])
        wtmp <- tmp_grad[2:(length(tmp_grad) - 1)]
        D[ind[2:(length(ind) - 1)]] <- dtmp * (cumsum(wtmp * xtmp) - xtmp * cumsum(wtmp) + xtmp * sum(wtmp * (1 - xtmp)))
    }
}

# generate output
res <- list("psi" = psi, "L" = L, "conc" = conc, "D" = D)
return(res)

}
