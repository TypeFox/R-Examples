LocalExtend <- function(X, isKnot, X2, psi2){

# extrapolates psi2 linearly between knots

KK <- which(isKnot == 1)
psi <- rep(0, length(X))
psi[KK] <- psi2

for (k in 1:(length(KK) - 1)){
    if (KK[k + 1] > (KK[k] + 1)){
        ind <- (KK[k] + 1):(KK[k + 1] - 1)
        lambda <- (X[ind] - X2[k]) / (X2[k + 1] - X2[k])

        # psi("not knot") = linearly interpolated from neighboring knots
        psi[ind] <- (1 - lambda) * psi2[k] + lambda * psi2[k + 1]
    }
}

return(psi)
}
