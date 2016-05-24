survexp_plot <-
function(futime,status,age,sex,entry_date,ratetable=survexp.fr,
               main="Observed and expected survival",xlab="Time (years)",ylab="Survival",
               col.km="black",lwd.km=2,lty.km=1,conf.int.km=TRUE,
               col.exp="blue",lwd.exp=2,lty.exp=1,
               main.rel="Relative survival",ylab.rel="Relative survival",
               col.rel="black",lwd.rel=2,lty.rel=1,times=seq(0,max(futime,na.rm=TRUE)/365.241,length=6)[-1],alpha=0.05,
               xscale=365.241,...){
  data=na.omit(data.frame(futime,status,age,sex,entry_date))
  par(mfrow=c(1,2))
  # survie observée
  km=survfit(Surv(futime,status)~1,data=data)
  plot(km,main=main,xlab=xlab,ylab=ylab,lwd=lwd.km,col=col.km,conf.int=conf.int.km,xscale=xscale,...)  
  # survie attendue
  expected=survexp(futime~1,rmap=list(sex=sex,year=entry_date,age=age),conditional=TRUE,data=data,ratetable=ratetable)
  lines(expected,lwd=lwd.exp,col=col.exp,lty=lty.exp,xscale=xscale)
  # legend
  legend("bottomleft",col=c(col.km,col.exp),lty=c(lty.km,lty.exp),lwd=c(lwd.km,lwd.exp),legend=c("Observed","Expected"))
  # relative
  fit=AER(data$futime,data$status,data$age,data$sex,data$entry_date,PY.stand=1,ratetable=ratetable,alpha=alpha)
  coef=log(fit$AER)
  coef.lo=log(fit$AER.lo)
  coef.up=log(fit$AER.up)    
  surv.rel=exp(-exp(coef)*times)
  surv.rel.lo=exp(-exp(coef.lo)*times)
  surv.rel.up=exp(-exp(coef.up)*times)
  fun_curve=function(x){exp(-exp(coef)*x)}
  curve(fun_curve,from=0,to=1.1*max(times),main=main.rel,xlab=xlab,ylab=ylab.rel,
        col=col.rel,lwd=lwd.rel,lty=lty.rel,ylim=c(0,1))    
  arrows(times,surv.rel.lo,times,surv.rel.up,col=col.rel,lwd=lwd.rel,angle=90,code=3,length=0.07)
  points(times,surv.rel,pch=19,col=col.rel)
  return(cbind(times,surv.rel,surv.rel.lo,surv.rel.up))  
}
