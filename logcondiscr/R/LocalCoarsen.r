LocalCoarsen <- function(X, W, isKnot){

# Re-distributes X and W on knot points
#
# OUTPUT
# X2 : X(isKnot == 1)
# W2 : W re-weighted on W(isKnot == 1) 

KK <- which(isKnot == 1)
X2 <- X[KK]
W2 <- W[KK]
for (k in 1:(length(KK) - 1)){

    if (KK[k + 1] > (KK[k] + 1)){

        # indices of observations between two knots
        ind <- (KK[k] + 1) : (KK[k + 1] - 1)
        
        # relative distance of these observations to left knot
        lambda <- (X[ind] - X2[k]) / (X2[k + 1] - X2[k])

        # weights of "no knots" re-distributed to closest "is knot"
        W2[k]     <- W2[k]     + sum(W[ind] * (1 - lambda))
        W2[k + 1] <- W2[k + 1] + sum(W[ind] * lambda)
    } # end if
} # end for

res <- list("X2" = X2, "W2" = W2)
return(res)
}
