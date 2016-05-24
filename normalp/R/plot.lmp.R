plot.lmp<-function(x,...){
dat<-x
subc<-deparse(dat$call)
plot(dat$fitted,dat$resid,main="",xlab="Fitted values",ylab="Residuals")
abline(h=0,lty=3)
mtext("Residuals vs Fitted", 3, 0.25)
title(sub=subc)
op<-par(ask = dev.interactive())
qqnorm(x$rs,main="",ylab="Standardized residuals",xlab="Theoretical Quantiles")
mtext("Normal Q-Q plot", 3, 0.25)
title(sub=subc)
rs<-paramp(dat$residuals,p=dat$p)
res<-(dat$residuals-rs$mp)/rs$sp
qqnormp(res,p=rs$p,main="Exponential Power Distribution Q-Q plot",ylab="Standardized residuals",xlab="Theoretical Quantiles")
title(sub=subc)
plot(dat$fitted,(abs(res))^(1/rs$p),main="",xlab="Fitted values",ylab=expression(sqrt("Standardized residuals",p)))
mtext("Scale-Location plot", 3, 0.25)
title(sub=subc)
par(op)   
}

