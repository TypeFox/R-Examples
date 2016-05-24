library("freqdom")
library(MASS)
library(mvtnorm)

data(gsl)
n = dim(precip$coef)[2]
ntest = floor(n*0.2)
n = n-ntest
T = 100
freq = pi*-T:T/T
X = t(precip$coef)
#PC = prcomp(X)

Y = t(elev$coef)
#PC = prcomp(Y)

K = 4 #4 #speclagreg.K(X,Y,freq)

# Estimate operators
Afull = speclagreg(X[1:n,],Y[1:n,],lags=0,freq=freq,K=K)
#R = reglag.boot(X[1:n,],Y[1:n,],Afull,rep=5,plot=TRUE,K=K)
#A = timedom.trunc(Afull,-1:1)
#W = reglag.significance(X[1:n,],Y[1:n,],A, alpha = 0.05, plot=TRUE)

# Notice no significance - this indicates that
# the linear regression (lagged and nonlagged) is irrelevant in this case
A = Afull

Yest = A %c% X
print(paste("Relative error: ", MSE(Y[1:ntest + n,],Yest[1:ntest + n,]) / MSE(Y[1:ntest + n,],0)))

print(TsiegiM(X[1:ntest + n,],Y[1:ntest + n,],Afull,lags=NULL))

