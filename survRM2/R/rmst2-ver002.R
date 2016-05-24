######################################
# rmst2 sample data
######################################
rmst2.sample.data=function(){

data("pbc", envir = environment()) ;
tmp <- get("pbc", envir  = environment()) ;

D=tmp[1:312,c(2:5,10,11,13,19)]
D$time=D$time/365.25 
D$status=as.numeric(D$status==2);
D$arm=as.numeric(D$trt==1) ;
DA=D[,c(1, 2, 9, 4:8)]
DA
}

######################################
# rmst2reg (Lu) -- hidden 
######################################
rmst2reg=function(y, delta, x, arm, tau, type="difference", alpha=0.05)
{

if(type!="difference" && type!="ratio" && type!="lossratio")
  print("Type must be difference, ratio or lossratio.")

if(type=="difference" || type=="ratio" || type=="lossratio"){

n=length(y)
x=cbind(1, x)
p=length(x[1,])

y0=pmin(y, tau)
d0=delta
d0[y0==tau]=1

d10=d0[arm==1]
d00=d0[arm==0]
y10=y0[arm==1]
y00=y0[arm==0]
x1=x[arm==1,]
x0=x[arm==0,]
n1=length(d10)
n0=length(d00)


id1=order(y10)
y10=y10[id1]
d10=d10[id1]
x1=x1[id1,]

id0=order(y00)
y00=y00[id0]
d00=d00[id0]
x0=x0[id0,]

fitc1=survfit(Surv(y10, 1-d10)~1)
fitc0=survfit(Surv(y00, 1-d00)~1)

weights1=d10/rep(fitc1$surv, table(y10))
weights0=d00/rep(fitc0$surv, table(y00))

weights=c(weights1, weights0)

if(type=="difference")
  {fitt=lm(c(y10,y00)~rbind(x1, x0)-1, weights=weights)
   beta0=fitt$coef
   
   error1=y10-as.vector(x1%*%beta0)
   score1=x1*weights1*error1

   error0=y00-as.vector(x0%*%beta0)
   score0=x0*weights0*error0
   }

if(type=="ratio")
  {fitt=glm(c(y10,y00)~rbind(x1, x0)-1, family="quasipoisson", weights=weights)
   beta0=fitt$coef
   
   error1=y10-exp(as.vector(x1%*%beta0))
   score1=x1*weights1*error1

   error0=y00-exp(as.vector(x0%*%beta0))
   score0=x0*weights0*error0
   }


if(type=="lossratio")
  {fitt=glm(c(tau-y10,tau-y00)~rbind(x1, x0)-1, family="quasipoisson", weights=weights)
   beta0=fitt$coef
   
   error1=tau-y10-exp(as.vector(x1%*%beta0))
   score1=x1*weights1*error1

   error0=tau-y00-exp(as.vector(x0%*%beta0))
   score0=x0*weights0*error0
   }



kappa.arm1=matrix(0, n1, p)
for(i in 1:n1)
   {kappa1=score1[i,]

    kappa2=apply(score1[y10>=y10[i],,drop=F], 2, sum)*(1-d10[i])/sum(y10>=y10[i])

    kappa3=rep(0, p)
    for(k in 1:n1)
       { if(y10[k]<=y10[i])
           kappa3=kappa3+apply(score1[y10>=y10[k],,drop=F], 2, sum)*(1-d10[k])/(sum(y10>=y10[k]))^2
       }

    kappa.arm1[i,]=kappa1+kappa2-kappa3
    }


kappa.arm0=matrix(0, n0, p)
for(i in 1:n0)
   {kappa1=score0[i,]

    kappa2=apply(score0[y00>=y00[i],,drop=F], 2, sum)*(1-d00[i])/sum(y00>=y00[i])

    kappa3=rep(0, p)
    for(k in 1:n0)
       {if(y00[k]<=y00[i])
           kappa3=kappa3+apply(score0[y00>=y00[k],,drop=F], 2, sum)*(1-d00[k])/(sum(y00>=y00[k]))^2
       }

    kappa.arm0[i,]=kappa1+kappa2-kappa3
    }


if(type=="difference")
  {gamma=t(kappa.arm1)%*%kappa.arm1+t(kappa.arm0)%*%kappa.arm0
   A=t(x)%*%x
   varbeta=solve(A)%*%gamma%*%solve(A)
   }

if(type=="ratio" || type=="lossratio")
  {gamma=t(kappa.arm1)%*%kappa.arm1+t(kappa.arm0)%*%kappa.arm0
   A=t(x*exp(as.vector(x%*%beta0)))%*%x
   varbeta=solve(A)%*%gamma%*%solve(A)
   }  


if(type=="difference")
  {beta0=beta0
   se0=sqrt(diag(varbeta))
   z0=beta0/se0
   p0=1-pchisq(z0^2, 1)
   cilow=beta0-se0*qnorm(1-alpha/2)
   cihigh=beta0+se0*qnorm(1-alpha/2)
   result=cbind(beta0, se0, z0, p0, cilow, cihigh)
   colnames(result)=c("coef", "se(coef)", "z","p",paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""))
   }

if(type=="ratio" || type=="lossratio")
  {beta0=beta0
   se0=sqrt(diag(varbeta))
   z0=beta0/se0
   p0=1-pchisq(z0^2, 1)
   r0=exp(beta0)
   cilow=exp(beta0-se0*qnorm(1-alpha/2))
   cihigh=exp(beta0+se0*qnorm(1-alpha/2))
   result=cbind(beta0, se0, z0, p0, exp(beta0), cilow, cihigh)
   colnames(result)=c("coef", "se(coef)", "z","p","exp(coef)",paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""))
   }

if(p==2)
   rownames(result)=c("intercept", "x")

if(p>2)
   rownames(result)=c("intercept", colnames(x[,-1]))

return(result)
}
}


#############################
# rmst1 (one-arm) -- hidden 
#############################
rmst1=function(time, status, tau, alpha=0.05){
 #-- time
 #-- statuts
 #-- tau -- truncation time
 #-- alpha -- gives (1-alpha) confidence interval

  ft= survfit(Surv(time, status)~1)
  idx=ft$time<=tau

  wk.time=sort(c(ft$time[idx],tau))
  wk.surv=ft$surv[idx]
  wk.n.risk =ft$n.risk[idx]
  wk.n.event=ft$n.event[idx]

  time.diff <- diff(c(0, wk.time))
  areas <- time.diff * c(1, wk.surv)
  rmst = sum(areas)
  rmst

  wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0, 
		        wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
  wk.var =c(wk.var,0)
  rmst.var = sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
  rmst.se  = sqrt(rmst.var)
       
  #--- check ---
  # print(ft, rmean=tau)

  #--- output ---
  out=matrix(0,2,4)
  out[1,]=c(rmst, rmst.se, rmst-qnorm(1-alpha/2)*rmst.se, rmst+qnorm(1-alpha/2)*rmst.se)
  out[2,]=c(tau-out[1,1], rmst.se, tau-out[1,4], tau-out[1,3])
  rownames(out)=c("RMST","RMTL") 
  colnames(out)=c("Est.", "se", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""))
  
  Z=list()
  Z$result=out
  Z$rmst = out[1,]
  Z$rmtl = out[2,]
  Z$tau=tau
  Z$rmst.var = rmst.var
  Z$fit=ft
  class(Z)="rmst1"

  return(Z)  

}



#########################################
# rmst2 (2-arm) contrast (main function)
#########################################
#rmst2=function(time, status, arm, tau=NULL, covariates=NULL, adjust.method="reg", alpha=0.05){
 rmst2=function(time, status, arm, tau=NULL, covariates=NULL,                      alpha=0.05){
 #-- time
 #-- statuts
 #-- arm (1 or 0)
 #-- covariates (matrix)
 #-- adjust = "reg"-- regression ("aug" -- augumentation)
 #-- alpha=0.05

#==================================
#  initial check 
#==================================
	
#--- tau ---
idx=arm==0; tt=time[idx]; event=status[idx] ; tau0=max(tt[event==1]) ; tau0max=max(tt)
idx=arm==1; tt=time[idx]; event=status[idx] ; tau1=max(tt[event==1]) ; tau1max=max(tt)
tau_default = min(tau0,    tau1)
tau_max     = min(tau0max, tau1max)


#---------------------
if(!is.null(tau)){
	if(tau<= tau_default){
       NOTE=paste("The truncation time: tau =", tau, " was specified.")
       }
	if(tau>= tau_default & tau<= tau_max){
       NOTE=paste("The truncation time: tau =", tau, " was specified, but there are no observed events after tau=," ,tau, "on either or both groups. Make sure that the size of riskset at tau=," ,tau, "is large enough in each group.")
		}	
	if(                    tau>= tau_max){
     stop(paste("The truncation time, tau, needs to be shorter than or equal to the minimum of the largest observed time on each of the two groups: ", round(tau_max, digits=3)))
        }
}
#---------------------
if(is.null(tau)){
   tau = tau_default
   NOTE=(paste("The truncation time, tau, was not specified. Thus, the default tau (the minimum of the largest observed event time on each of the two groups)", round(tau_default, digits=3)," is used."))
}
#---------------------



Z=list()
Z$tau=tau
Z$note=NOTE

#==================================
#  unadjusted analysis 
#==================================
if(is.null(covariates)){

wk1=rmst1(time[arm==1], status[arm==1], tau)
wk0=rmst1(time[arm==0], status[arm==0], tau)

Z$RMST.arm1=wk1
Z$RMST.arm0=wk0

  
#--- contrast (RMST difference) ---
 rmst.diff.10     = wk1$rmst[1]-wk0$rmst[1]
 rmst.diff.10.se  = sqrt(wk1$rmst.var + wk0$rmst.var)
 rmst.diff.10.low = rmst.diff.10 - qnorm(1-alpha/2)*rmst.diff.10.se
 rmst.diff.10.upp = rmst.diff.10 + qnorm(1-alpha/2)*rmst.diff.10.se
 rmst.diff.pval   = pnorm(-abs(rmst.diff.10)/rmst.diff.10.se)*2 
 rmst.diff.result = c(rmst.diff.10, rmst.diff.10.low, rmst.diff.10.upp, rmst.diff.pval)

#--- contrast (RMST ratio) ---
 rmst.log.ratio.10     = log(wk1$rmst[1]) - log(wk0$rmst[1])
 rmst.log.ratio.10.se  = sqrt(wk1$rmst.var/wk1$rmst[1]/wk1$rmst[1] + wk0$rmst.var/wk0$rmst[1]/wk0$rmst[1])
 rmst.log.ratio.10.low = rmst.log.ratio.10 - qnorm(1-alpha/2)*rmst.log.ratio.10.se
 rmst.log.ratio.10.upp = rmst.log.ratio.10 + qnorm(1-alpha/2)*rmst.log.ratio.10.se
 rmst.log.ratio.pval   = pnorm(-abs(rmst.log.ratio.10)/rmst.log.ratio.10.se)*2 
 rmst.ratio.result     = c(exp(rmst.log.ratio.10), exp(rmst.log.ratio.10.low), exp(rmst.log.ratio.10.upp),rmst.log.ratio.pval)

#--- contrast (RMTL ratio  0/1) ---
# rmtl.log.ratio.01     = log(wk0$rmtl[1]) - log(wk1$rmtl[1])
# rmtl.log.ratio.01.se  = sqrt(wk1$rmst.var/wk1$rmtl[1]/wk1$rmtl[1] + wk0$rmst.var/wk0$rmtl[1]/wk0$rmtl[1])
# rmtl.log.ratio.01.low = rmtl.log.ratio.01 - qnorm(1-alpha/2)*rmtl.log.ratio.01.se
# rmtl.log.ratio.01.upp = rmtl.log.ratio.01 + qnorm(1-alpha/2)*rmtl.log.ratio.01.se
# rmtl.log.ratio.pval   = pnorm(-abs(rmtl.log.ratio.01)/rmtl.log.ratio.01.se)*2 
# rmtl.ratio.result     = c(exp(rmtl.log.ratio.01), exp(rmtl.log.ratio.01.low), exp(rmtl.log.ratio.01.upp),rmtl.log.ratio.pval)

#--- contrast (RMTL ratio  1/0) ---
 rmtl.log.ratio.10     = log(wk1$rmtl[1]) - log(wk0$rmtl[1])
 rmtl.log.ratio.10.se  = sqrt(wk1$rmst.var/wk1$rmtl[1]/wk1$rmtl[1] + wk0$rmst.var/wk0$rmtl[1]/wk0$rmtl[1])
 rmtl.log.ratio.10.low = rmtl.log.ratio.10 - qnorm(1-alpha/2)*rmtl.log.ratio.10.se
 rmtl.log.ratio.10.upp = rmtl.log.ratio.10 + qnorm(1-alpha/2)*rmtl.log.ratio.10.se
 rmtl.log.ratio.pval   = pnorm(-abs(rmtl.log.ratio.10)/rmtl.log.ratio.10.se)*2 
 rmtl.ratio.result     = c(exp(rmtl.log.ratio.10), exp(rmtl.log.ratio.10.low), exp(rmtl.log.ratio.10.upp),rmtl.log.ratio.pval)



 #--- results ---
 out=rbind(rmst.diff.result, rmst.ratio.result , rmtl.ratio.result )
 rownames(out)=c("RMST (arm=1)-(arm=0)","RMST (arm=1)/(arm=0)","RMTL (arm=1)/(arm=0)") 
 colnames(out)=c("Est.", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""), "p")
  

 #--- output ---
 Z$unadjusted.result = out
 # Z$RMST.difference=out[1,]
 # Z$RMST.ratio=out[2,]
 # Z$RMTL.ratio=out[3,]
 # Z
}
	
#==================================
#  Adjusted analysis 
#==================================
if (!is.null(covariates)){

## if (adjust.method=="reg"){

aa=rmst2reg(time, status, as.matrix(cbind(arm, covariates)), arm, tau, type="difference", alpha=alpha)
bb=rmst2reg(time, status, as.matrix(cbind(arm, covariates)), arm, tau, type="ratio", alpha=alpha)
cc=rmst2reg(time, status, as.matrix(cbind(arm, covariates)), arm, tau, type="lossratio", alpha=alpha)


 #--- output ---
 out.adj=rbind(aa[2,c(1,5,6,4)], bb[2, c(5,6,7,4)], cc[2, c(5,6,7,4)])
 rownames(out.adj)=c("RMST (arm=1)-(arm=0)","RMST (arm=1)/(arm=0)","RMTL (arm=1)/(arm=0)") 
 colnames(out.adj)=c("Est.", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""), "p")
  

 #--- output ---
 Z$adjusted.result = out.adj

 Z$RMST.difference.adjusted = aa
 Z$RMST.ratio.adjusted      = bb
 Z$RMTL.ratio.adjusted      = cc

## }else{
##  stop "Please sepcify adjust.method" 	
## }

}


 class(Z)="rmst2"
 
 Z

}


######################################
# print.rmst2 (hidden)
######################################
print.rmst2=function(x, digits=3, ...){

cat("\n")

cat(x$note,"\n\n")
   


#--- unadjusted analysis --
 if(!is.null(x$unadjusted.result)){

 	RMST=rbind(x$RMST.arm1$result[1,], x$RMST.arm0$result[1,])
 	RMTL=rbind(x$RMST.arm1$result[2,], x$RMST.arm0$result[2,])
    rownames(RMST)=c("RMST (arm=1)","RMST (arm=0)")
    rownames(RMTL)=c("RMTL (arm=1)","RMTL (arm=0)")

	cat ("Restricted Mean Survival Time (RMST) by arm \n")

    prmatrix(round(RMST , digits=digits))

    cat("\n\n")

	cat ("Restricted Mean Time Lost (RMTL) by arm \n")

    prmatrix(round(RMTL, digits=digits))

    cat("\n\n")

	cat ("Between-group contrast \n")

  	prmatrix(round(x$unadjusted.result, digits=digits))

 }


#--- Adjusted analysis --
 if(!is.null(x$adjusted.result)){

	cat ("Summary of between-group contrast (adjusted for the covariates) \n")
  	prmatrix(round(x$adjusted.result, digits=digits))

    cat("\n\n")

	cat ("Model summary (difference of RMST) \n")
 	prmatrix(round(x$RMST.difference.adjusted, digits=digits))

    cat("\n\n")

	cat ("Model summary (ratio of RMST) \n")
 	prmatrix(round(x$RMST.ratio.adjusted, digits=digits))

    cat("\n\n")

	cat ("Model summary (ratio of time-lost) \n")
 	prmatrix(round(x$RMTL.ratio.adjusted, digits=digits))


 }

 invisible(x)
}


######################################
# plot.rmst2 (hidden)
######################################
plot.rmst2=function(x, xlab="", ylab="", col="red", col.RMST="pink", col.RMTL="orange",density=80, angle=85,...){

 if(is.null(x$unadjusted.result)) stop("Please run rmst2 without covariates")
 
 if(!is.null(x$unadjusted.result)){
  
  ft1=x$RMST.arm1$fit
  ft0=x$RMST.arm0$fit
  tau=x$tau
 
  par(mfrow=c(1,2))

  #=== arm 1 ===
  fit=ft1
 
  tmp.xx=c(0, fit$time); tmp.yy=c(1, fit$surv) ;
  idx=tmp.xx<=tau
  y.tau = min(tmp.yy[idx])
  xx=c(tmp.xx[idx],   tau)
  yy=c(tmp.yy[idx], y.tau)  
  x.step=sort(c(0, xx, xx))
  y.step=rev(sort(c(1,1,yy, yy[-length(yy)])))

  #--------
  plot(fit, mark.time=F, conf.int=F, lwd=2, main="arm=1", xlab=xlab, ylab=ylab, col=col, sub=paste("RMST:",round(x$RMST.arm1$rmst[1], digits=2)))

  for (i in seq(1, length(x.step), by=2)){  
  polygon(c(x.step[i], x.step[i+1], x.step[i+1], x.step[i]), c(0, 0, y.step[i+1], y.step[i]), col= col.RMST, density=density, angle=angle, lwd=2)
  }
  for (i in seq(1, length(x.step), by=2)){  
  polygon(c(x.step[i], x.step[i+1], x.step[i+1], x.step[i]), c(y.step[i], y.step[i+1], 1,1), col= col.RMTL, density=density, angle=angle, lwd=2)
  }

  x.step=sort(c(0, tmp.xx, tmp.xx))
  y.step=rev(sort(c(1,1,tmp.yy, tmp.yy[-length(tmp.yy)])))
  lines(x.step, y.step, col=col, lwd=3) 
  # text(5,0.4, paste(round(rmst$rmst[1], digits=2),"years"), cex=0.9)


  #=== arm 0 ===
  fit=ft0
 
  tmp.xx=c(0, fit$time); tmp.yy=c(1, fit$surv) ;
  idx=tmp.xx<=tau
  y.tau = min(tmp.yy[idx])
  xx=c(tmp.xx[idx],   tau)
  yy=c(tmp.yy[idx], y.tau)  
  x.step=sort(c(0, xx, xx))
  y.step=rev(sort(c(1,1,yy, yy[-length(yy)])))
 
  #--------
  plot(fit, mark.time=F, conf.int=F, lwd=2, main="arm=0", xlab=xlab, ylab=ylab, col=col, sub=paste("RMST:",round(x$RMST.arm0$rmst[1], digits=2)))
  for (i in seq(1, length(x.step), by=2)){  
  polygon(c(x.step[i], x.step[i+1], x.step[i+1], x.step[i]), c(0, 0, y.step[i+1], y.step[i]), col= col.RMST, density=density, angle=angle, lwd=2)
  }
  for (i in seq(1, length(x.step), by=2)){  
  polygon(c(x.step[i], x.step[i+1], x.step[i+1], x.step[i]), c(y.step[i], y.step[i+1], 1,1), col= col.RMTL, density=density, angle=angle, lwd=2)
  }

  x.step=sort(c(0, tmp.xx, tmp.xx))
  y.step=rev(sort(c(1,1,tmp.yy, tmp.yy[-length(tmp.yy)])))
  lines(x.step, y.step, col=col, lwd=3) 
  # text(5,0.4, paste(round(rmst$rmst[1], digits=2),"years"), cex=0.9)

 }


}


######################################
# Example (PBC)
######################################
if(0){
D=rmst2.sample.data()
time=D$time
status=D$status
arm=D$arm
tau=NULL
x=D[,c(4,6,7)]

rmst2(time, status, arm)
rmst2(time, status, arm, tau=12)

a=rmst2(time, status, arm, tau=10)
print(a)
plot(a, xlab="Years", ylab="Probability", density=60)

a=rmst2(time, status, arm, tau=10, covariates=x)
print(a)
a=rmst2(time, status, arm, tau=10, covariates=x, alpha=0.1)
print(a)
}