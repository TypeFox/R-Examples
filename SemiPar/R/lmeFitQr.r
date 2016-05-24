###### R-function: lmeFitQr ##########

# Computes the QR decomposition for a
# lme() fit with pre-specified covariance
# structure.

# Last changed: 04 AUG 2004

lmeFitQr <- function(y,X,Z,G,resid.var)
{
   p <- ncol(X)
   q <- ncol(Z)

   # Form the "C matrix" and then obtain R^(-1/2)C

   C.mat <- cbind(X,Z)

   R.neghalf.C <- C.mat/sqrt(resid.var)

   # Augment R^{-1/2} C with zeros and G^{-1/2} to get C.a

   C.augment <- matrix(0,p+q,p+q)
   
   G.neghalf <- backsolve(chol(G),diag(rep(1,q)))

   C.augment[(p+1):(p+q),(p+1):(p+q)] <- G.neghalf

   C.a <- rbind(R.neghalf.C,C.augment)

   y.a <- c(y,rep(0,p+q))/sqrt(resid.var)

   # compute conditional covariance matrix

   Ca.qr <- lm.fit(C.a,y.a)

   return(Ca.qr)
}

########## End of lmeFitQr ##########
