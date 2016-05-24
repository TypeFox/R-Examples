## ----echo=FALSE----------------------------------------------------------
set.seed(22381950)

## ------------------------------------------------------------------------
n <- 1000
Xlogreg <- data.frame(x1=rbinom(n,1,0.4), x2=runif(n,0,1), x3=rnorm(n))
xbeta <- -1+Xlogreg$x1+Xlogreg$x2+Xlogreg$x3
Xlogreg$y <- rbinom(n,1,exp(xbeta)/(1+exp(xbeta)))

## ------------------------------------------------------------------------
library(spcadjust)
chartlogreg <- new("SPCCUSUM",model=SPCModellogregLikRatio(Delta= 1, formula="y~x1+x2+x3"))
xihat <- xiofdata(chartlogreg,Xlogreg)
xihat

## ------------------------------------------------------------------------
cal <- SPCproperty(data=Xlogreg,
            nrep=100,chart=chartlogreg,
            property="calARL",params=list(target=1000),quiet=TRUE)
cal

## ----echo=FALSE----------------------------------------------------------
set.seed(2238195)

## ------------------------------------------------------------------------
n <- 100
newXlogreg <- data.frame(x1=rbinom(n,1,0.4), x2=runif(n,0,1), x3=rnorm(n))
newxbeta <- -1+newXlogreg$x1+newXlogreg$x2+newXlogreg$x3
newXlogreg$y <- rbinom(n,1,exp(newxbeta)/(1+exp(newxbeta)))
S <- runchart(chartlogreg, newdata=newXlogreg,xi=xihat)

## ----fig=TRUE,fig.width=10,fig.height=4,echo=FALSE-----------------------
par(mfrow=c(1,1),mar=c(4,5,0,0))
plot(S,ylab=expression(S[t]),xlab="t",type="b",ylim=range(S,cal@res+1,cal@raw))
lines(c(0,100),rep(cal@res,2),col="red")
lines(c(0,100),rep(cal@raw,2),col="blue")
legend("topleft",c("Adjusted Threshold","Unadjusted Threshold"),col=c("red","blue"),lty=1)

## ----echo=FALSE----------------------------------------------------------
set.seed(22)

## ------------------------------------------------------------------------
n <- 100
newXlogreg <- data.frame(x1=rbinom(n,1,0.4), x2=runif(n,0,1), x3=rnorm(n))
outind <- c(rep(0,50),rep(1,50))
newxbeta <- -1+newXlogreg$x1+newXlogreg$x2+newXlogreg$x3+outind
newXlogreg$y <- rbinom(n,1,exp(newxbeta)/(1+exp(newxbeta)))
S <- runchart(chartlogreg, newdata=newXlogreg,xi=xihat)

## ----fig=TRUE,fig.width=10,fig.height=4,echo=FALSE-----------------------
par(mfrow=c(1,1),mar=c(4,5,0,0))
plot(S,ylab=expression(S[t]),xlab="t",type="b",ylim=range(S,cal@res+1,cal@raw))
lines(c(0,100),rep(cal@res,2),col="red")
lines(c(0,100),rep(cal@raw,2),col="blue")
legend("topleft",c("Adjusted Threshold","Unadjusted Threshold"),col=c("red","blue"),lty=1)

