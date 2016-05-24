## ----echo=FALSE----------------------------------------------------------
set.seed(12381900)

## ------------------------------------------------------------------------
X <-  rnorm(250)

## ----fig=TRUE,fig.width=8,fig.height=3,echo=FALSE------------------------
par(mar=c(4,5,0,0))
plot(-(250:1),X,xlab="t",ylab=expression(X[t]))

## ------------------------------------------------------------------------
library(spcadjust)
chart <- new("SPCCUSUM",model=SPCModelNormal(Delta=1));
xihat <- xiofdata(chart,X)
str(xihat)

## ------------------------------------------------------------------------
cal <- SPCproperty(data=X,nrep=50,
            property="calARL",chart=chart,params=list(target=100),quiet=TRUE)
cal

## ------------------------------------------------------------------------
newX <- rnorm(100)
S <- runchart(chart, newdata=newX,xi=xihat)

## ----fig=TRUE,fig.width=10,fig.height=4----------------------------------
par(mfrow=c(1,2),mar=c(4,5,0,0))
plot(newX,xlab="t")
plot(S,ylab=expression(S[t]),xlab="t",type="b",ylim=range(S,cal@res+1,cal@raw))
lines(c(0,100),rep(cal@res,2),col="red")
lines(c(0,100),rep(cal@raw,2),col="blue")
legend("topleft",c("Adjusted Threshold","Unadjusted Threshold"),col=c("red","blue"),lty=1)

## ------------------------------------------------------------------------
newX <- rnorm(100,mean=c(rep(0,50),rep(1,50)))
S <- runchart(chart, newdata=newX,xi=xihat)

## ----fig=TRUE,fig.width=10,fig.height=4,echo=FALSE-----------------------
par(mfrow=c(1,2),mar=c(4,4,0,0))
plot(newX,xlab="t")
plot(S,ylab=expression(S[t]),xlab="t",type="b",ylim=pmin(range(S,cal@res,cal@raw),15))
lines(c(0,100),rep(cal@res,2),col="red")
lines(c(0,100),rep(cal@raw,2),col="blue")
legend("topleft",c("Adjusted Threshold","Unadjusted Threshold"),col=c("red","blue"),lty=1)

