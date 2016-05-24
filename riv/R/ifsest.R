# 11) Function to get the sample IF of an S-est.

# Input:
# Z is a (nxr) data matrix.
# (L,V) is the location/scatter S-estimator for this data matrix
# R is the number of MonteCarlo replications to compute expectations
# bdp=Breakdown point to compute S-estimator and c-constant of
# Tukey's biweight function
# c=constant for Tukey's biweight function
# Econst=list containing expectations needed for the IF of the S-est.
IF.Sest <- function(Z, L, V, c, Econst) {
    B <- t(chol(V))
    Binv <- solve(B)

    z0 <- Binv %*% (t(Z) - L)
    norm.z0 <- sqrt(colSums(z0^2))
    
    IF.M0 <- t(z0) * (1/Econst$nu) * psi.bisquare(norm.z0, c)
    IFL <- B %*% t(IF.M0)
    
    IFV <- apply(z0, 2, IF.obs, B, Econst, c)
    dim(IFV) <- c(nrow(z0), nrow(z0), ncol(z0))
    
    list(IFL=IFL, IFV=IFV)   
}
