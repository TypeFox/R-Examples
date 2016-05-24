## ------------------------------------------------------------------------
library(spcadjust)
model <- SPCModelNormal(Delta=1)
model$Pofdata
model$Pofdata <- function(data){
      list(mu= median(data), sd= mad(data), m=length(data))
}

## ------------------------------------------------------------------------
X <-  rnorm(100)
chartrobust <- new("SPCCUSUM",model=model)
SPCproperty(data=X,nrep=50,property="calARL",
            chart=chartrobust,params=list(target=100),quiet=TRUE)

## ------------------------------------------------------------------------
SPCModelExponential=function(Delta=1){
    structure(
        list(
            Pofdata=function(data){
                list(lambda=1/mean(data), n=length(data))
            },
            xiofP=function(P) P$lambda,
            resample=function(P) rexp(P$n,rate=P$lambda),
            getcdfupdates=function(P, xi) {
                function(x){ if(Delta<1)
                                 pmax(0,1-exp(-P$lambda*(x-log(Delta))/(xi*(1-Delta))))
                else
                    pmin(1,exp(-P$lambda*(log(Delta)-x)/(xi*(Delta-1))))
                         }
            },
            updates=function(xi,data) log(Delta)-xi*(Delta-1)*data),
        class="SPCDataModel")
}

ExpCUSUMchart=new("SPCCUSUM",model=SPCModelExponential(Delta=1.25))

## ----echo=FALSE----------------------------------------------------------
set.seed(223819)

## ----fig=TRUE,fig.width=8,fig.height=4-----------------------------------
X <- rexp(1000)
plot(runchart(ExpCUSUMchart, newdata=rexp(100),xi=xiofdata(ExpCUSUMchart,X)),
     ylab=expression(S[t]),xlab="t",type="b")

## ------------------------------------------------------------------------
SPCproperty(data=X,nrep=100,property="hitprob",
chart=ExpCUSUMchart,params=list(threshold=1,nsteps=100),
covprob=c(0.5,0.9),quiet=TRUE)
SPCproperty(data=X,nrep=100,property="ARL",
            chart=ExpCUSUMchart,params=list(threshold=3),covprob=c(0.5,0.9),quiet=TRUE)
SPCproperty(data=X,nrep=100,property="calARL",chart=ExpCUSUMchart,
            params=list(target=1000),covprob=c(0.5,0.9),quiet=TRUE) 

## ----fig=TRUE,fig.width=10,fig.height=4----------------------------------
cal <- SPCproperty(data=X,nrep=1000,property="calARL",chart=ExpCUSUMchart,
                   params=list(target=1000),quiet=TRUE,parallel=1)
S <- runchart(ExpCUSUMchart, newdata=rexp(100),xi=xiofdata(ExpCUSUMchart,X))
par(mfrow=c(1,1),mar=c(4,5,0,0))
plot(S,ylab=expression(S[t]),xlab="t",type="b",ylim=range(S,cal@res+1,cal@raw))
lines(c(0,100),rep(cal@res,2),col="red")
lines(c(0,100),rep(cal@raw,2),col="blue")
legend("topleft",c("Adjusted Threshold","Unadjusted Threshold"),col=c("red","blue"),lty=1)

