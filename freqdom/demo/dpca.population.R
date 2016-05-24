library('freqdom')

if (!requireNamespace("MARSS", quietly = TRUE)) {
  stop("MARSS package is needed for this demo to work. Please install it.",
    call. = FALSE)
}

library('MARSS')

data('kestrel')

# prepare data
X = kestrel[,-1] 
m = colMeans(X)
X = t(t(X) - m)
n = dim(X)[1]

## Static PCA ##
PR = prcomp(X)
Y1 = PR$x
Y1[,-1] = 0
Xpca = Y1 %*% t(PR$rotation)

## Dynamic PCA ##
XI.est = dprcomp(X,q=20,weights="Daniell")  # finds the optimal filter
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

# PLOT 10 observations reconstructed from the first component
ind = 1:10 + 20
par(mfrow=c(1,3))
ylim=c(-1,1)
matplot(X[ind,],ylim=ylim,xlab="Intraday time", ylab="Sqrt(PM10)",t='l')
matplot(Xpca[ind,],ylim=ylim,xlab="Intraday time", ylab="Sqrt(PM10)",t='l')
matplot(Xdpca.est[ind,],ylim=ylim,xlab="Intraday time", ylab="Sqrt(PM10)",t='l')
par(mfrow=c(1,1))
