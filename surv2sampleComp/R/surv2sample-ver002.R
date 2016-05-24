#############################
# surv2sampleComp (ver002) add average t-years and percentiles
#############################
q95=function(x){ quantile(x, prob=c(0.025, 0.975))}

surv2sample=function(time, status, arm, npert=1000, timepoints=c(12,24,36,40), quanprobs =c(0.1, 0.15, 0.2), tau=NULL, SEED=NULL, procedure="KM"){
	
if(!is.null(SEED)){ set.seed(SEED) }

if(is.null(tau)){
  idx=arm==0; tt=time[idx]; event=status[idx] ; tau0=max(tt[event==1])
  idx=arm==1; tt=time[idx]; event=status[idx] ; tau1=max(tt[event==1])
  tau=min(tau0,tau1)
}

if(procedure=="KM"){
  idx=arm==0; km0=KM2.pert(time[idx], status[idx], npert=npert, timepoints=timepoints, quanprobs = quanprobs, tau=tau)
  idx=arm==1; km1=KM2.pert(time[idx], status[idx], npert=npert, timepoints=timepoints, quanprobs = quanprobs, tau=tau)
}
if(procedure=="GG"){
  idx=arm==0; km0=GG2.boot(time[idx], status[idx], npert=npert, timepoints=timepoints, quanprobs = quanprobs, tau=tau)
  idx=arm==1; km1=GG2.boot(time[idx], status[idx], npert=npert, timepoints=timepoints, quanprobs = quanprobs, tau=tau)
}
  

#--- each arm ---
# wk0=cbind(km0$rmst, tau-km0$rmst, km0$tyearprobs, km0$percentiles)
# wk1=cbind(km1$rmst, tau-km1$rmst, km1$tyearprobs, km1$percentiles)
  wk0=cbind(km0$rmst, tau-km0$rmst, km0$tyearprobs, km0$percentiles, km0$tyearprobs.ave, km0$percentiles.ave)
  wk1=cbind(km1$rmst, tau-km1$rmst, km1$tyearprobs, km1$percentiles, km1$tyearprobs.ave, km1$percentiles.ave)
  se0=apply(wk0[-1,], 2, sd)
  se1=apply(wk1[-1,], 2, sd)
  ci0=apply(wk0[-1,], 2, q95)
  ci1=apply(wk1[-1,], 2, q95)

# measure=c("RMST","Loss time", paste("Prob at",timepoints), paste("Quantile at", quanprobs*100,"%")) 
  measure=c("RMST","Loss time", paste("Prob at",timepoints), paste("Quantile at", quanprobs*100,"%"), "Ave of t-year event rates","Ave percentiles") 
  
  out.group0=cbind(t(wk0[1,]), t(ci0), se0)
  out.group1=cbind(t(wk1[1,]), t(ci1), se1)
  rownames(out.group0)=measure
  rownames(out.group1)=measure
  colnames(out.group0)=c("Est.","Lower 95%","Upper 95%","SE")
  colnames(out.group1)=c("Est.","Lower 95%","Upper 95%","SE")
  
#--- contrast ---
 K=ncol(wk0)
 outq=c()
 for (k in 1:K){
  q0=wk0[,k] ; q1=wk1[,k]
  diff01=q0-q1 
  diff10=q1-q0 
  ratio01=q0/q1
  ratio10=q1/q0
  outq=cbind(outq, diff01, diff10, ratio01, ratio10)
 }

 measures=rep(measure, each=4) 
 contrast=rep(c("Group0-Group1","Group1-Group0","Group0/Group1","Group1/Group0"), K)

#======================================
if(procedure=="KM"){


 #--- p-val --
 se=apply(outq[-1,], 2, sd)
 pval=pnorm(-abs(outq[1,])/se)*2 
 
 idx=contrast=="Group0/Group1"|contrast=="Group1/Group0"
 for (j in 1:length(pval)){
   if(contrast[j]=="Group0/Group1" | contrast[j]=="Group1/Group0"){
     for (m in 1:length(measure)){
     	if(measures[j]==measure[m]){
           se[j]=sqrt( (se0[m]/out.group0[m,1])^2 + (se1[m]/out.group1[m,1])^2)
        }
      }
    }
 }
 pval[idx]=pnorm(-abs(log(outq[1,idx]))/se[idx])*2 

 #--- confidence interals --
 lower=outq[1,] - 1.96*se
 upper=outq[1,] + 1.96*se
 idx=contrast=="Group0/Group1"|contrast=="Group1/Group0"
 lower[idx]=exp(log(outq[1,idx]) - 1.96*se[idx])
 upper[idx]=exp(log(outq[1,idx]) + 1.96*se[idx])



 #--- results ---
 out.contrast=cbind(outq[1,], lower, upper, pval)
 rownames(out.contrast)=paste(measures, contrast)
 colnames(out.contrast)=c("Est.","Lower 95%","Upper 95%", "p-val")
 inf.method="Perturbation resampling"

}
#======================================
if(procedure=="GG"){

 #--- confidence interals (bootstrap percentile) --
 cband=apply(outq[-1,], 2, q95)
 lower=cband[1,]
 upper=cband[2,]

 out.contrast=cbind(outq[1,], lower, upper)
 rownames(out.contrast)=paste(measures, contrast)
 colnames(out.contrast)=c("Est.","Lower 95%","Upper 95%")
 inf.method="Bootstrap percentile method"

}
#======================================


 #--- output ---
 Z=list()
 Z$procedure = procedure
 Z$method = inf.method
 Z$survfit = survfit(Surv(time, status)~arm)
 Z$tau=tau
 Z$npert=npert 
 Z$timepoints=timepoints
 Z$quanprobs=quanprobs
 Z$contrast.all=out.contrast
 Z$group0=out.group0
 Z$group1=out.group1
 Z$RMST=out.contrast[measures=="RMST",]
 Z$RMLT=out.contrast[measures=="Loss time",]
 Z$contrast.diff10=out.contrast[contrast=="Group1-Group0",]
 Z$contrast.diff01=out.contrast[contrast=="Group0-Group1",]
 Z$contrast.ratio01=out.contrast[contrast=="Group0/Group1",]
 Z$contrast.ratio10=out.contrast[contrast=="Group1/Group0",]
 
 class(Z)="surv2sample"
 return(Z)

}


######################################
# Example (E4A03) (confidential)
######################################
# mydata<-read.dta("../../CHECK/data-e4a03-covs.dta")
# time=mydata$time
# status=mydata$status
# arm=mydata$arm
# comp=surv2sample(time, status, arm, npert=200, timepoints=c(12,24,36,40), quanprobs =c(0.05, 0.1, 0.15, 0.2), SEED=1223, tau=40)
# comp$RMST
# comp$contrast.ratio10
# comp$contrast.ratio01
# attributes(comp)
# print(comp)
