library("freqdom")
library("mvtnorm")
data(pm10)

n = dim(X$coef)[2]

Xar = t(X$coef[,1:n]) 
Yar = t(Y$coef[,1:n])

p = prcomp(Xar,center=TRUE,retx=TRUE)
Xar = p$x[,1:5]
p = prcomp(Yar,center=TRUE,retx=TRUE)
Yar = p$x[,1:5]

ntest = floor(n * 0.1)
n = n - ntest
T = 100
freq = pi*-T:T/T

K = speclagreg.K(Xar[1:n,],(Yar[1:n,]),freq)

#plot.alpha(t(X$coef),1:4)

# Estimate the regressors and check significance
Afull = speclagreg(Xar[1:n,],(Yar[1:n,]),lags=-1:1,K=K)
#A = lagreg.est(Xar[1:n,],matrix(Yar[1:n,]),lags=0:1,K=10)
#R = reglag.boot(Xar[1:n,],(Xar[1:n,]),Afull,rep=1,plot=FALSE,K=K)
W = reglag.significance(Xar[1:n,],(Yar[1:n,]),Afull, alpha = 0.05, plot=FALSE)

# Lag zero is significant. Let's check errors in practice. First all lags
A = Afull
Yest = (A %c% Xar)
Ytest = matrix(Yar[1:ntest + n,])
Ytest.est = matrix(Yest[1:ntest + n,])
print(paste("Relative error: ", MSE(Ytest,Ytest.est) / MSE(Ytest,0)))

print(TsiegiM(Xar[1:ntest + n,],Yar[1:ntest + n,],Afull,lags=NULL))
