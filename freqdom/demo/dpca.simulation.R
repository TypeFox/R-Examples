library("freqdom")

set.seed(1)
n = 100

# We illustrate here an extreme case of very strong
# temporal dependance in the process
# Generate n observations of a AR(1) process
P = matrix(c(1,-9,-3,1),2,2)
P = 0.9* P / norm.spec(P)
X = rar(n, Psi=P, noise=function(n){0.1*rnorm(n)})

# standard PCA
PR = prcomp(X)
PR$x = X %*% t(PR$rotation)
PR$x[,2] = 0 # use only the first component
X.est.static = PR$x %*% PR$rotation

# dynamic PCA
XI = dprcomp(X,-10:10,q=15)
Y = XI %c% X
Y[,2] = 0 # use only the first score
X.est = t(rev(XI)) %c% Y # recover X from first score

# Compute variance explained by PCA
M1 = MSE(X,X.est.static)
M = MSE(X,0)

cat("NMSE PCA =  ")
cat(M1/M) # approximation error of the first score

# Compute variance explained by DPCA
M1 = MSE(X,X.est)
cat("\nNMSE DPCA = ")
cat(M1/M) # approximation error of the first score
cat("\n")

# Plot scores
par(mfrow=c(1,3))

ind = 1:20

matplot(X[ind,],t='l',xlab="time")
title("Original time series")
matplot(X.est.static[ind,],xlab="time",t='l')
title("PCA recovered time series")
matplot(X.est[ind,],xlab="time",t='l')
title("DPCA recovered time series")
par(mfrow=c(1,1))
