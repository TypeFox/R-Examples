if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("MASS package is needed for this demo to work. Please install it.",
    call. = FALSE)
}

library("freqdom")
library(MASS)
library(mvtnorm)

data(ice.river)
n = dim(ice.river)[1]
ntest = floor(n*0.2)
n = n-ntest
T = 100
freq = pi*-T:T/T
rnames = rownames(ice.river)
X = as.matrix(ice.river[,3:4])
MX = colMeans(X)
X = t(t(X) - MX)
#PC = prcomp(X)

Y = as.matrix(ice.river[,1:2])
MY = colMeans(Y)
Y = t(t(Y) - MY)
#PC = prcomp(Y)

K = 4 #4 #speclagreg.K(X,Y,freq)

# Estimate operators
Afull = speclagreg(X[1:n,],Y[1:n,],lags=-10:10,freq=freq,K=K)
#R = reglag.boot(X[1:n,],Y[1:n,],Afull,rep=5,plot=TRUE,K=K)
#A = timedom.trunc(Afull,-1:1)
#W = reglag.significance(X[1:n,],Y[1:n,],A, alpha = 0.05, plot=TRUE)

# Notice no significance - this indicates that
# the linear regression (lagged and nonlagged) is irrelevant in this case
A = Afull

Yest = (A %c% X)
cat(paste("Relative error: ", MSE(Y[1:ntest + n,],Yest[1:ntest + n,]) / MSE(Y[1:ntest + n,],0)),"\n")

Y = t(t(Y) + MY)
Yest = t(t(Yest) + MY)
plot(Yest[,2],t='l',ylim = c(0,150))
lines(Y[,2],t='l',col='red')
