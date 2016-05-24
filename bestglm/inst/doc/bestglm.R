### R code from vignette source 'bestglm.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
online <- FALSE ## if set to FALSE the local copy of MSFT.rda
                ## is used instead of get.hist.quote()
options(prompt = "R> ")


###################################################
### code chunk number 2: bestglm-zprostate
###################################################
library(bestglm)
data(zprostate)


###################################################
### code chunk number 3: PRINTzprostate
###################################################
str(zprostate)


###################################################
### code chunk number 4: bestglm-train
###################################################
train<-(zprostate[zprostate[,10],])[,-10]
X<-train[,1:8]
y<-train[,9]
out <- summary(regsubsets(x = X, y = y, nvmax=ncol(X)))
Subsets <- out$which
RSS <- out$rss


###################################################
### code chunk number 5: print1
###################################################
cbind(as.data.frame(Subsets), RSS=RSS)


###################################################
### code chunk number 6: bestglmbicqArgs
###################################################
args(bestglm)


###################################################
### code chunk number 7: bestglmbicqOut
###################################################
Xy<-cbind(as.data.frame(X), lpsa=y)
out<-bestglm(Xy)
names(out)


###################################################
### code chunk number 8: bestglm-aic
###################################################
bestglm(Xy, IC="AIC")


###################################################
### code chunk number 9: bestglm-bic
###################################################
bestglm(Xy, IC="BIC")


###################################################
### code chunk number 10: bestglm-bicg
###################################################
bestglm(Xy, IC="BICg")


###################################################
### code chunk number 11: bestglmNumericalExample
###################################################
set.seed(1233211235)
p<-5   #number of inputs
n<-100  #number of observations
X<-matrix(rnorm(n*p), ncol=p)
err<-rnorm(n)
y<- 0.1*(X[,1]+X[,2]+X[,3])+err
Xy<-as.data.frame(cbind(X,y))
names(Xy)<-c(paste("X",1:p,sep=""),"y")
ans <- bestglm(Xy)
ans$Subsets


###################################################
### code chunk number 12: bestglmNumericalExampleContinued
###################################################
ans$qTable


###################################################
### code chunk number 13: bestglmbicq
###################################################
data(zprostate)
train<-(zprostate[zprostate[,10],])[,-10]
X<-train[,1:8]
y<-train[,9]
Xy<-cbind(as.data.frame(X), lpsa=y)
out <- bestglm(Xy, IC="BICq")


###################################################
### code chunk number 14: bestglmDeleted
###################################################
set.seed(123321123)
bestglm(Xy, IC="CV", t=10)


###################################################
### code chunk number 15: bestglmfolds
###################################################
set.seed(2377723)
ind<-sample(rep(1:10,length=67))
ind


###################################################
### code chunk number 16: bestglmfoldsPiOne
###################################################
(1:67)[1==ind]


###################################################
### code chunk number 17: bestglmfolds
###################################################
tabulate(ind)


###################################################
### code chunk number 18: bestglmCVHTF
###################################################
set.seed(2377723)
out<-bestglm(Xy, IC="CV", CVArgs=list(Method="HTF", K=10, REP=1)) 
out


###################################################
### code chunk number 19: plotCVHTF (eval = FALSE)
###################################################
## cverrs <- out$Subsets[,"CV"]
## sdCV<- out$Subsets[,"sdCV"]
## CVLo <- cverrs - sdCV
## CVHi<- cverrs + sdCV
## ymax <- max(CVHi)
## ymin <- min(CVLo)
## k <- 0:(length(cverrs)-1)
## plot(k, cverrs, xlab="Subset Size", ylab="CV Error", ylim=c(ymin,ymax), type="n", yaxt="n")
## points(k, cverrs, cex=2, col="red", pch=16)
## lines(k, cverrs, col="red", lwd=2)
## axis(2, yaxp=c(0.6,1.8,6))
## segments(k, CVLo, k, CVHi, col="blue", lwd=2)
## eps <- 0.15
## segments(k-eps, CVLo, k+eps, CVLo, col="blue", lwd=2)
## segments(k-eps, CVHi, k+eps, CVHi, col="blue", lwd=2)
## #cf. oneSDRule
## indBest <- oneSdRule(out$Subsets[,c("CV", "sdCV")])
## abline(v=indBest-1, lty=2)
## indMin <- which.min(cverrs)
## fmin <- sdCV[indMin]
## cutOff <- fmin + cverrs[indMin]
## abline(h=cutOff, lty=2)
## indMin <- which.min(cverrs)
## fmin <- sdCV[indMin]
## cutOff <- fmin + cverrs[indMin]
## min(which(cverrs < cutOff))


###################################################
### code chunk number 20: plotCVHTF-repeat
###################################################
cverrs <- out$Subsets[,"CV"]
sdCV<- out$Subsets[,"sdCV"]
CVLo <- cverrs - sdCV
CVHi<- cverrs + sdCV
ymax <- max(CVHi)
ymin <- min(CVLo)
k <- 0:(length(cverrs)-1)
plot(k, cverrs, xlab="Subset Size", ylab="CV Error", ylim=c(ymin,ymax), type="n", yaxt="n")
points(k, cverrs, cex=2, col="red", pch=16)
lines(k, cverrs, col="red", lwd=2)
axis(2, yaxp=c(0.6,1.8,6))
segments(k, CVLo, k, CVHi, col="blue", lwd=2)
eps <- 0.15
segments(k-eps, CVLo, k+eps, CVLo, col="blue", lwd=2)
segments(k-eps, CVHi, k+eps, CVHi, col="blue", lwd=2)
#cf. oneSDRule
indBest <- oneSdRule(out$Subsets[,c("CV", "sdCV")])
abline(v=indBest-1, lty=2)
indMin <- which.min(cverrs)
fmin <- sdCV[indMin]
cutOff <- fmin + cverrs[indMin]
abline(h=cutOff, lty=2)
indMin <- which.min(cverrs)
fmin <- sdCV[indMin]
cutOff <- fmin + cverrs[indMin]
min(which(cverrs < cutOff))


###################################################
### code chunk number 21: bestglmCVDH
###################################################
set.seed(2377723)
bestglm(Xy, IC="CV", CVArgs=list(Method="DH", K=10, REP=1))
bestglm(Xy, IC="CV", CVArgs=list(Method="DH", K=10, REP=1)) 
bestglm(Xy, IC="CV", CVArgs=list(Method="DH", K=10, REP=1))  


###################################################
### code chunk number 22: bestglmLOOCV
###################################################
bestglm(Xy, IC="LOOCV")


###################################################
### code chunk number 23: bestglmManpowerAIC
###################################################
data(manpower)
bestglm(manpower, IC="AIC")


###################################################
### code chunk number 24: bestglmManpowerBIC
###################################################
bestglm(manpower, IC="BIC")


###################################################
### code chunk number 25: bestglmManpowerBICg
###################################################
bestglm(manpower, IC="BICg")


###################################################
### code chunk number 26: bestglmManpowerBICghalf
###################################################
bestglm(manpower, IC="BICg", t=0.5)


###################################################
### code chunk number 27: bestglmManpowerBICq
###################################################
out<-bestglm(manpower, IC="BICq")
out


###################################################
### code chunk number 28: bestglmManpowerBestq
###################################################
out$Bestq


###################################################
### code chunk number 29: bestglmManpowerSubsets
###################################################
out$Subsets


###################################################
### code chunk number 30: bestglmManpowerqTable
###################################################
out$qTable


###################################################
### code chunk number 31: bestglmSAheartBinomialFULL
###################################################
data(SAheart)
out<-bestglm(SAheart, IC="BICq", t=1, family=binomial)
out


###################################################
### code chunk number 32: bestglmSAheartBinomialStep
###################################################
ans<-glm(chd~., data=SAheart)
q<-0.25
n<-nrow(SAheart)
k<-log(n) - 2*log(q/(1-q))
step(ans, k=k)


###################################################
### code chunk number 33: bestglmSAheartBinomialGaussian
###################################################
out<-bestglm(SAheart, IC="BICq", t=0.25)
out$Subsets


###################################################
### code chunk number 34: bestglmNuclear
###################################################
data(znuclear)
bestglm(znuclear, IC="AIC")


###################################################
### code chunk number 35: bestglmDetroit1
###################################################
data(Detroit)
X<-as.data.frame(scale(Detroit[,c(1,2,4,6,7,11)]))
y<-Detroit[,ncol(Detroit)]
Xy<-cbind(X,HOM=y)
out <- lm(HOM~., data=Xy)
step(out, k=log(nrow(Xy)))


###################################################
### code chunk number 36: bestglmDetroit2
###################################################
out<-bestglm(Xy, IC="BIC")
out


###################################################
### code chunk number 37: bestglmDetroit2
###################################################
out$qTable


###################################################
### code chunk number 38: bestglmDetroit3
###################################################
bestglm(Xy,IC="BICq", t=0.05)


###################################################
### code chunk number 39: bestglmDetroit3
###################################################
bestglm(Xy,IC="BICq", t=0.00005)


###################################################
### code chunk number 40: bestglmDetroit4
###################################################
set.seed(1233211)
bestglm(Xy, IC="CV", CVArgs=list(Method="d", K=4, REP=50))


###################################################
### code chunk number 41: bestglmAirQualityFullModel
###################################################
data(AirQuality)
bestglm(AirQuality,IC="BICq",t=1)


###################################################
### code chunk number 42: bestglmAirQualityAIC
###################################################
bestglm(AirQuality,IC="AIC")


###################################################
### code chunk number 43: bestglmFiresAIC
###################################################
data(Fires)
bestglm(Fires, IC="AIC")


###################################################
### code chunk number 44: bestglmManpowerNullRegression
###################################################
set.seed(123312123)
X<-as.data.frame(matrix(rnorm(50), ncol=2, nrow=25))
y<-rnorm(25)
Xy<-cbind(X, y=y)
bestglm(Xy)


###################################################
### code chunk number 45: bestglmLogistic
###################################################
set.seed(231231)
n<-500
K<-10 #number of inputs not counting constant
a<--1
b<-c(c(9,6,4,8)/3, rep(0, K-4))
X<-matrix(rnorm(n*K), ncol=K)
L<-a+X%*%b
p<-plogis(L)
Y<-rbinom(n=n, size=1, prob=p)
X<-as.data.frame(X)
#X<-as.matrix.data.frame(X)
out<-glm(Y~., data=X, family=binomial)
summary(out)


###################################################
### code chunk number 46: bestglmBinomial
###################################################
set.seed(231231)
n<-500
K<-8 #number of inputs not counting constant
m<-100 #binomial - number of trials
a<-2
b<-c(c(9,6,4,8)/10, rep(0, K-4))
X<-matrix(rnorm(n*K), ncol=K)
L<-a+X%*%b
p<-plogis(L)
Y<-rbinom(n=n, size=m, prob=p)
Y<-cbind(Y, m-Y)
dimnames(Y)[[2]]<-c("S","F")
X<-as.data.frame(X)
out<-glm(Y~., data=X, family=binomial)
summary(out)


###################################################
### code chunk number 47: bestglmBinomialBestModel
###################################################
Xy<-cbind(X, Y)
bestglm(Xy, family=binomial)


###################################################
### code chunk number 48: bestglmBinomialFactorVariable
###################################################
set.seed(33344111)
n<-500
K<-4 #number of quantitative inputs not counting constant
m<-100 #binomial - number of trials
a<-2 #intercept
dayNames<-c("Sunday","Monday","Tuesday","Wednesday","Friday","Saturday")
Days<-data.frame(d=factor(rep(dayNames, n))[1:n])
Xdays<-model.matrix(~d, data=Days)
bdays<-c(7,2,-7,0,2,7)/10
Ldays<-Xdays%*%bdays
b<-c(c(9,6)/10, rep(0, K-2))
X<-matrix(rnorm(n*K), ncol=K)
L<-a+X%*%b
L<-L + Ldays
p<-plogis(L)
Y<-rbinom(n=n, size=m, prob=p)
Y<-cbind(Y, m-Y)
dimnames(Y)[[2]]<-c("S","F")
X<-as.data.frame(X)
X<-data.frame(X, days=Days)
out<-glm(Y~., data=X, family=binomial)
anova(out,test="Chisq")


###################################################
### code chunk number 49: bestglmBinomialFactorVariableBestModel
###################################################
Xy <- cbind(X, Y)
out<-bestglm(Xy, IC="BICq", family=binomial)
out


###################################################
### code chunk number 50: bestglmPoissonSimulate
###################################################
set.seed(231231)
n<-500
K<-4 #number of inputs not counting constant
a<- -1
b<-c(c(1,0.5), rep(0, K-2))
X<-matrix(rnorm(n*K), ncol=K)
L<-a+X%*%b
lambda<-exp(L)
Y<-rpois(n=n, lambda=lambda)
X<-as.data.frame(X)
#X<-as.matrix.data.frame(X)
out<-glm(Y~., data=X, family=poisson)
summary(out)


###################################################
### code chunk number 51: bestglmPoissonFit
###################################################
Xy <- data.frame(X, y=Y)
bestglm(Xy, family=poisson)


###################################################
### code chunk number 52: bestglmGammaSimulateAndFit
###################################################
GetGammaParameters<-function(muz, sdz){
    phi<-(sdz/muz)^2
    nu<-1/phi
    lambda<-muz/nu
    list(shape=nu, scale=lambda)
    }
set.seed(321123)
test<-rnorm(20)
n<-500
b<-c(0.25, 0.5, 0, 0)
b0<-0.3
K<-length(b)
sdz<-1
X<-matrix(rnorm(n*K), ncol=K)
L<-b0+X%*%b
muHat<-exp(L)
gp <- GetGammaParameters(muHat, sdz)
zsim<-rgamma(n, shape=gp$shape, scale=gp$scale)
Xy<-data.frame(as.data.frame.matrix(X), y=zsim)
out<-glm(y~., data=Xy, family=Gamma(link=log))
summary(out)


###################################################
### code chunk number 53: bestglmGammaBest
###################################################
bestglm(Xy, family=Gamma(link=log))


###################################################
### code chunk number 54: bestglmShao
###################################################
`testCorrect` <-
function(ans,NB){
NBfit<-names(coef(ans))[-1]
ans<-ifelse(length(NBfit)==length(NB)&(!any(is.na(match(NBfit,NB)))),1,0)
ans
}
#
NSIM<-5 #100 simulations takes about 12 sec
data(Shao)
set.seed(123321123)
X<-as.matrix.data.frame(Shao)
BETA<-list(b1=c(0,0,4,0),b2=c(0,0,4,8),b3=c(9,0,4,8),b4=c(9,6,4,8))
NamesBeta<-list(b1=c("x4"), b2=c("x4","x5"), b3=c("x2","x4","x5"),b4=c("x2","x3","x4","x5"))
hitsBIC<-hitsEBIC<-hitsQBIC<-numeric(4)
startTime<-proc.time()[1]
for (iB in 1:4){
    b<-BETA[[iB]]
    NB<-NamesBeta[[iB]]
    for (iSIM in 1:NSIM){
        y <- 2+X%*%b+rnorm(40)
        Xy<-cbind(Shao, y)
        hitsBIC[iB]<-hitsBIC[iB]+testCorrect(bestglm(Xy, IC="BIC")$BestModel,NB)
        hitsEBIC[iB]<-hitsEBIC[iB]+testCorrect(bestglm(Xy, IC="BICg")$BestModel,NB)
        hitsQBIC[iB]<-hitsQBIC[iB]+testCorrect(bestglm(Xy, IC="BICq")$BestModel,NB)
    }
}
endTime<-proc.time()[1]
totalTime<-endTime-startTime
ans<-matrix(c(hitsBIC,hitsEBIC,hitsQBIC),byrow=TRUE,ncol=4)
dimnames(ans)<-list(c("BIC","BICg","BICq"),1:4)
ans<-t(ans)/NSIM
ans
totalTime



###################################################
### code chunk number 55: bestglmSIM
###################################################
set.seed(123321123) 
startTime<-proc.time()[1]
#NSIM<-50
NSIM<-5
p<-25  #number of inputs
n<-30  #number of observations
ans<-numeric(4)
names(ans)<-c("AIC", "BIC", "BICg", "BICq")
for (iSim in 1:NSIM){
    X<-matrix(rnorm(n*p), ncol=p)
    y<-rnorm(n)
    Xy<-as.data.frame(cbind(X,y))
    names(Xy)<-c(paste("X",1:p,sep=""),"y")
    bestAIC <- bestglm(Xy, IC="AIC")
    bestBIC <- bestglm(Xy, IC="BIC")
    bestEBIC <- bestglm(Xy, IC="BICg")
    bestQBIC <- bestglm(Xy, IC="BICq", t=0.05)
    ans[1] <- ans[1] +  length(coef(bestAIC$BestModel))-1
    ans[2] <- ans[2] +  length(coef(bestBIC$BestModel))-1
    ans[3] <- ans[3] +  length(coef(bestEBIC$BestModel))-1
    ans[4] <- ans[4] +  length(coef(bestQBIC$BestModel))-1
}
endTime<-proc.time()[1]
totalTime<-endTime-startTime
totalTime
ans


