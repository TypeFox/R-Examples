library(sme)
data(MTB)
fits <- sme(MTB)
sapply(fits,logLik)

fits <- sme(MTB,numerOfThreads=10)
sapply(fits,logLik)

fit <- sme(MTB[MTB$variable==6031,c("y","tme","ind")])
logLik(fit)

data(inflammatory)
fit <- sme(inflammatory,deltaNM=0.1,deltaEM=0.1)
logLik(fit)

fit <- sme(inflammatory,knots=c(29.5,57,84.5),deltaEM=0.1,deltaNM=0.1)
AIC(fit)
AICc(fit)
BIC(fit)
BICn(fit)