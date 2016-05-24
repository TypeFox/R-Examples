## ----echo=FALSE----------------------------------------------------------
set.seed(12381900)

## ------------------------------------------------------------------------
n <- 1000
Xlinreg <- data.frame(x1= rbinom(n,1,0.4), x2= runif(n,0,1), x3= rnorm(n))
Xlinreg$y <- 2 + Xlinreg$x1 + Xlinreg$x2 + Xlinreg$x3 + rnorm(n)

## ------------------------------------------------------------------------
library(spcadjust)
chartlinreg <- new("SPCCUSUM",model=SPCModellm(Delta=1,formula="y~x1+x2+x3"))
xihat <- xiofdata(chartlinreg,Xlinreg)
xihat

## ------------------------------------------------------------------------
cal <- SPCproperty(data=Xlinreg,
             nrep=100,
             property="calARL",chart=chartlinreg,params=list(target=100),quiet=TRUE)
cal

## ----echo=FALSE----------------------------------------------------------
set.seed(12381903)

## ------------------------------------------------------------------------
n <- 100
newXlinreg <- data.frame(x1= rbinom(n,1,0.4), x2= runif(n,0,1), x3=rnorm(n))
newXlinreg$y <- 2 + newXlinreg$x1 + newXlinreg$x2 + newXlinreg$x3 + rnorm(n)
S <- runchart(chartlinreg, newdata=newXlinreg,xi=xihat)

## ----fig=TRUE,fig.width=10,fig.height=4,echo=FALSE-----------------------
par(mfrow=c(1,1),mar=c(4,5,0,0))
plot(S,ylab=expression(S[t]),xlab="t",type="b",ylim=range(S,cal@res+1,cal@raw))
lines(c(0,100),rep(cal@res,2),col="red")
lines(c(0,100),rep(cal@raw,2),col="blue")
legend("topleft",c("Adjusted Threshold","Unadjusted Threshold"),col=c("red","blue"),lty=1)

## ------------------------------------------------------------------------
n <- 100
newXlinreg <- data.frame(x1= rbinom(n,1,0.4), x2= runif(n,0,1),x3=rnorm(n))
outind <- c(rep(0,50),rep(1,50))
newXlinreg$y <- 2 + newXlinreg$x1 + newXlinreg$x2 + newXlinreg$x3 + rnorm(n)+outind
S <- runchart(chartlinreg, newdata=newXlinreg,xi=xihat)

## ----fig=TRUE,fig.width=10,fig.height=4,echo=FALSE-----------------------
par(mfrow=c(1,1),mar=c(4,5,0,0))
plot(S,ylab=expression(S[t]),xlab="t",type="b",ylim=range(S,cal@res+1,cal@raw))
lines(c(0,100),rep(cal@res,2),col="red")
lines(c(0,100),rep(cal@raw,2),col="blue")
legend("topleft",c("Adjusted Threshold","Unadjusted Threshold"),col=c("red","blue"),lty=1)

