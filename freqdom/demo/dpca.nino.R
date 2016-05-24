library("freqdom")

if (!requireNamespace("fda", quietly = TRUE)) {
  stop("fda package is needed for this demo to work. Please install it.",
    call. = FALSE)
}
if (!requireNamespace("tseries", quietly = TRUE)) {
  stop("tseries package is needed for this demo to work. Please install it.",
    call. = FALSE)
}

library("tseries")
library("fda")

data(nino)

args=seq(0,1,length=12)
#basis = create.fourier.basis(rangeval=c(0, 1), nbasis=9)
basis = create.bspline.basis(rangeval=c(0, 1), nbasis=9, norder=3)

nino3 = matrix(nino3[1:(12*48)],12)
nino3.4 = matrix(nino3.4[1:(12*48)],12)
X.fd = center.fd(Data2fd(args,nino3,basis))
Y.fd = center.fd(Data2fd(args,nino3.4,basis))
X = t(X.fd$coef)
Y = t(Y.fd$coef)

# Estimate operators
A = speclagreg(X,Y,lags=-5:5)

n = dim(X)[1]
## Static PCA ##
PR = prcomp(X)
Y1 = PR$x
Y1[,-1] = 0
Xpca = Y1 %*% t(PR$rotation)

## Dynamic PCA ##
V = inprod(basis,basis) # we need to take into account correlations between elements
XI.est = dprcomp(X,V=V,q=20,weights="Daniell")  # finds the optimal filter
Y.est = XI.est %c% X  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xdpca.est = t(rev(XI.est)) %c% Y.est   # deconvolution

# Write down results
ind = 1:n
cat("NMSE DPCA = ")
cat(MSE(X[ind,],Xdpca.est[ind,]) / MSE(X[ind,],0))
cat("\nNMSE PCA =  ")
cat(MSE(X[ind,],Xpca[ind,]) / MSE(X[ind,],0))
cat("\n")

# Creates functional objects
Xdpca.est.fd = fd(t(Re(Xdpca.est)),basis=X.fd$basis)
Xpca.fd = fd(t(Xpca),basis=X.fd$basis)

# PLOT 10 observations reconstructed from the first component
ind = 1:10 + 10
par(mfrow=c(1,3))
ylim = c(-1.5,1.5)
plot(X.fd[ind],ylim=ylim,xlab="Intraday time", ylab="Sqrt(PM10)")
title("Original curves")
plot(Xpca.fd[ind],ylim=ylim,xlab="Intraday time", ylab="Sqrt(PM10)")
title("PCA curves")
plot(Xdpca.est.fd[ind],ylim=ylim,xlab="Intraday time", ylab="Sqrt(PM10)")
title("DPCA curves")
par(mfrow=c(1,1))
