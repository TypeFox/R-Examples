roc.curves.boot <-
function(data,nb=1000,alfa=0.05,name,mod1,mod2,related=TRUE){  
  nomeg=paste(name)
  
  sim1.ind=unlist(data[1])
  sim2.ind=unlist(data[2])
  
  sim1.sta=unlist(data[3])
  sim2.sta=unlist(data[4])
  
  imax=length(sim1.sta[sim1.sta==0])
  jmax=length(sim1.sta[sim1.sta==1])
  imax2=length(sim2.sta[sim2.sta==0])
  jmax2=length(sim2.sta[sim2.sta==1])
  
  total_cases=imax+jmax
  total_cases2=imax2+jmax2
  sim1=t(rbind(sim1.sta,sim1.ind)[,order(sim1.sta,sim1.ind)])
  sim2=t(rbind(sim2.sta,sim2.ind)[,order(sim2.sta,sim2.ind)])
  dados=rbind(cbind(sim1[,2],sim1[,1]),cbind(sim2[,2],sim2[,1]))
  
  mod=2
  
  if(imax>imax2) maxImax = imax else maxImax = imax2
  if(jmax>jmax2) maxJmax = jmax else maxJmax = jmax2
  # norm defines the negative cases
  norm=array(data=NA,c(maxImax,mod))
  # abnorm defines the positive cases
  abnorm=array(data=NA,c(maxJmax,mod))
  
  # data are atributed to the arrays norm and abnorm
  for (i in 1:imax) {
    norm[i,1]=dados[i,1]
  }
  for (j in 1:jmax) {
    abnorm[j,1]=dados[j+imax,1]
  }
  for (i in 1:imax2) {
    norm[i,2]=dados[i+total_cases,1]
  }
  for (j in 1:jmax2) {
    abnorm[j,2]=dados[j+total_cases+imax2,1]
  }
  
  # calculates various values (frequencies, ...)
  sim1.pred <- prediction(sim1.ind, sim1.sta)
  sim2.pred <- prediction(sim2.ind, sim2.sta)
  
  # calculate points of the curves (TPR and FPR)
  sim1.curve <- performance(sim1.pred,"tpr","fpr")
  sim2.curve <- performance(sim2.pred,"tpr","fpr")
  
  # calculates the areas
  sim1.auc <- unlist(performance(sim1.pred,"auc")@y.values)
  sim2.auc <- unlist(performance(sim2.pred,"auc")@y.values)
  
  
  # calculates FPR and TPR
  sim1.fpr <- array(unlist(sim1.curve@x.values))
  sim1.tpr <- array(unlist(sim1.curve@y.values))
  sim2.fpr <- array(unlist(sim2.curve@x.values))
  sim2.tpr <- array(unlist(sim2.curve@y.values))
  
  diff.areas =sim1.auc-sim2.auc
  
  roc.curves.plot(sim1.curve,sim2.curve,mod1,mod2)
  
  #creation of the structure
  sim1.ind.sample=array(data=NA,c(length(sim1.sta),nb))
  sim2.ind.sample=array(data=NA,c(length(sim2.sta),nb))
  sim1.fpr.sample=array(data=NA,c(length(sim1.fpr),nb))
  sim2.fpr.sample=array(data=NA,c(length(sim2.fpr),nb))
  sim1.tpr.sample=array(data=NA,c(length(sim1.tpr),nb))
  sim2.tpr.sample=array(data=NA,c(length(sim2.tpr),nb))
  sim1.auc.sample=c()
  sim2.auc.sample=c()
  diff.auc.sample=c()
  dist1.sample=array(data=NA,c(100,nb))
  dist2.sample=array(data=NA,c(100,nb))
  diff.areas.sample=array(data=NA,c(101,nb))
  
  for(b in 1:nb){
    
    sim1.ind.sample[,b]=c(sample(norm[,1],imax,TRUE),sample(abnorm[,1],jmax,TRUE))
    sim2.ind.sample[,b]=c(sample(norm[,2],imax2,TRUE),sample(abnorm[,2],jmax2,TRUE))
    
    # calculates various values (frequencies, ...)
    sim1.pred.sample <- prediction(sim1.ind.sample[,b], sim1.sta)
    sim2.pred.sample <- prediction(sim2.ind.sample[,b], sim2.sta)
    
    # calculate points of the curves (TPR and FPR)
    sim1.curve.sample <- performance(sim1.pred.sample,"tpr","fpr")
    sim2.curve.sample <- performance(sim2.pred.sample,"tpr","fpr")
    
    # calculates the areas
    sim1.auc.sample[b] <- unlist(performance(sim1.pred.sample,"auc")@y.values)
    sim2.auc.sample[b] <- unlist(performance(sim2.pred.sample,"auc")@y.values)
    
    diff.auc.sample[b] <- sim1.auc.sample[b]-sim2.auc.sample[b]
    
    # calculate FPR and TPR
    sim1.fpr.sample[1:length(unlist(sim1.curve.sample@x.values)),b] <- array(unlist(sim1.curve.sample@x.values))
    sim1.tpr.sample[1:length(unlist(sim1.curve.sample@y.values)),b] <- array(unlist(sim1.curve.sample@y.values))
    sim2.fpr.sample[1:length(unlist(sim2.curve.sample@x.values)),b] <- array(unlist(sim2.curve.sample@x.values))
    sim2.tpr.sample[1:length(unlist(sim2.curve.sample@y.values)),b] <- array(unlist(sim2.curve.sample@y.values))
    
    result=rocsampling(array(unlist(sim1.curve.sample@x.values)),array(unlist(sim1.curve.sample@y.values)), array(unlist(sim2.curve.sample@x.values)),array(unlist(sim2.curve.sample@y.values)))
    dist1.sample[,b]=result$dist1
    dist2.sample[,b]=result$dist2
    diff.areas.sample[,b]=result$diffareas
  }
  
  
  result=rocsampling(sim1.fpr,sim1.tpr,sim2.fpr,sim2.tpr)
  rocsampling.summary(result,mod1,mod2)
  
  dist1=result$dist1                                
  dist2=result$dist2
  diff.dist=dist1-dist2
  diff.areas=result$diffareas
  
  rstar1=sqrt(total_cases)*(dist1.sample-dist1)
  rstar2=sqrt(total_cases)*(dist2.sample-dist2)
  rstar=rstar1-rstar2
  
  mean1.boot=c()
  sd1.boot=c()
  mean2.boot=c()
  sd2.boot=c()
  mean.boot=c()
  sd.boot=c()
  qL.boot=c()
  qU.boot=c()
  mean.rstar=c()
  sd.rstar=c()
  qL.rstar=c()
  qU.rstar=c()
  
  for(i in 1:length(dist1.sample[,2]))
  {
    mean1.boot[i]=mean(dist1.sample[i,])
    sd1.boot[i]=sd(dist1.sample[i,])
    mean2.boot[i]=mean(dist2.sample[i,])
    sd2.boot[i]=sd(dist2.sample[i,])
    mean.boot[i]=mean(dist1.sample[i,]-dist2.sample[i,])
    sd.boot[i]=sd(dist1.sample[i,]-dist2.sample[i,])
    qL.boot[i]=quantile(dist1.sample[i,]-dist2.sample[i,],alfa/2)
    qU.boot[i]=quantile(dist1.sample[i,]-dist2.sample[i,],1-alfa/2)
    mean.rstar[i]=mean(rstar[i,])
    sd.rstar[i]=sd(rstar[i,])
    qL.rstar[i]=quantile(rstar[i,],alfa/2)
    qU.rstar[i]=quantile(rstar[i,],1-alfa/2)
  }
  
  qL.t=qt(alfa/2,total_cases-2)
  qU.t=qt(1-alfa/2,total_cases-2)
  
  IC.par=array(data=NA,c(100,2))
  IC.per=array(data=NA,c(100,2))
  IC.bcp=array(data=NA,c(100,2))
  IC.par.rstar=array(data=NA,c(100,2))
  IC.per.rstar=array(data=NA,c(100,2))
  
  for(i in 1:length(qL.rstar))
  {
    IC.par[i,1]=(diff.dist[i]-mean.boot[i])+qL.t*sd.boot[i]
    IC.par[i,2]=(diff.dist[i]-mean.boot[i])+qU.t*sd.boot[i]
    IC.per[i,1]=(diff.dist[i]-mean.boot[i])+qL.boot[i]#diff.dist[i]+qL.boot[i]#*sd.boot[i]
    IC.per[i,2]=(diff.dist[i]-mean.boot[i])+qU.boot[i]#diff.dist[i]+qU.boot[i]#*sd.boot[i]
    IC.bcp[i,1]=(diff.dist[i]-mean.boot[i])+qL.boot[i]#diff.dist[i]+qL.boot[i]#*sd.boot[i]
    IC.bcp[i,2]=(diff.dist[i]-mean.boot[i])+qU.boot[i]#diff.dist[i]+qU.boot[i]#*sd.boot[i]
    IC.par.rstar[i,1]=diff.dist[i]+qL.t*sd.rstar[i]
    IC.par.rstar[i,2]=diff.dist[i]+qU.t*sd.rstar[i]
    IC.per.rstar[i,1]=diff.dist[i]+qL.rstar[i]#*sd.boot[i]
    IC.per.rstar[i,2]=diff.dist[i]+qU.rstar[i]#*sd.boot[i]
  }
  
  
  LB=min(c(IC.par[,1],IC.per[,1]))
  UB=max(c(IC.par[,2],IC.per[,2]))
  
  
  mean.areas.boot=c()
  sd.areas.boot=c()
  qL.areas.boot=c()
  qU.areas.boot=c()
  for(i in 1:length(diff.areas.sample[,2]))
  {
    mean.areas.boot[i]=mean(diff.areas.sample[i,])
    sd.areas.boot[i]=sd(diff.areas.sample[i,])            
    qL.areas.boot[i]=quantile(diff.areas.sample[i,],alfa/2)
    qU.areas.boot[i]=quantile(diff.areas.sample[i,],1-alfa/2)
  }
  IC.areas.par=array(data=NA,c(101,2))
  IC.areas.per=array(data=NA,c(101,2))
  for(i in 1:length(qL.areas.boot))
  {
    IC.areas.par[i,1]=(diff.areas[i]-mean.areas.boot[i])+qL.t*sd.areas.boot[i]
    IC.areas.par[i,2]=(diff.areas[i]-mean.areas.boot[i])+qU.t*sd.areas.boot[i]
    IC.areas.per[i,1]=(diff.areas[i]-mean.areas.boot[i])+qL.areas.boot[i]#diff.areas[i]+qL.areas.boot[i]#*sd.boot[i]
    IC.areas.per[i,2]=(diff.areas[i]-mean.areas.boot[i])+qU.areas.boot[i]#diff.areas[i]+qU.areas.boot[i]#*sd.boot[i]
  }
  
  LB=min(c(IC.areas.par[,1],IC.areas.per[,1]))
  UB=max(c(IC.areas.par[,2],IC.areas.per[,2]))
  
  degrees=result$lineslope*90/(pi/2)
  
  plot(degrees[1:101], c(rep(LB,50),rep(UB,51)), type='n', xlab="Degrees", ylab="Area between curves")
  abline(h=0, col = "gray60")
  lines(degrees[1:101], IC.areas.per[,1], col='red', lty=8)    
  #lines(degrees[1:101], IC.areas.par[,1], col='red', lty=1)    
  lines(degrees[1:101], diff.areas, col='blue')     
  lines(degrees[1:101], IC.areas.per[,2], col='green', lty=9)
  #lines(degrees[1:101], IC.areas.par[,2], col='green', lty=3)
  title("Areas Between ROC Curves")
  
  IC1.auc.per=c()
  IC2.auc.per=c()
  IC1.auc.per[1]=quantile(sim1.auc.sample,alfa/2)
  IC1.auc.per[2]=quantile(sim1.auc.sample,1-alfa/2)
  IC2.auc.per[1]=quantile(sim2.auc.sample,alfa/2)
  IC2.auc.per[2]=quantile(sim2.auc.sample,1-alfa/2)
  
  IC.diff.auc.per=c()
  IC.diff.auc.per[1]=quantile(diff.auc.sample,alfa/2)
  IC.diff.auc.per[2]=quantile(diff.auc.sample,1-alfa/2)
  
  resultdelong=comp.roc.delong(sim1.ind,sim1.sta,sim2.ind,sim2.sta,related)
  
  nc=0
  for(i in 1:(length(result$diffareas)-1))
  {
    if (result$diffareas[i]<0 && result$diffareas[i+1]>0) nc=nc+1
    if (result$diffareas[i]>0 && result$diffareas[i+1]<0) nc=nc+1
  }
  
  resultboot=comp.roc.curves(result,ci.flag=TRUE,graph.flag=T,nomeg)
  
  #if (k==1)
  #{
    #cat("\n   Delong & & & & & & & & Permutation & & & & Bootstrap & & & & & &\\\n")
    #cat("\n  Delong & & & & & & & & Permutation & & & & Bootstrap & & & & & &\\\n",file=paste(nomeg,"Results.txt"),append=TRUE)
    #cat("\n & AUC_1 & SE_1 & AUC_2 & SE_2 & R & Diff & Z & pvalue & AUC_1 & AUC_2 & pvalue & ncoss & IC_AUC_1 & IC_AUC_1 & IC_AUC_2 & IC_AUC_2 & IC_DIFF & IC_DIFF &\\\n")
    #cat("\n & AUC_1 & SE_1 & AUC_2 & SE_2 & R & Diff & Z & pvalue & AUC_1 & AUC_2 & pvalue & ncoss & IC_AUC_1 & IC_AUC_1 & IC_AUC_2 & IC_AUC_2 & IC_DIFF & IC_DIFF &\\\n",file=paste(nomeg,"Results_zang.txt"),append=TRUE);
  #}
  #cat(resultdelong$AUC[1],resultdelong$SE[1],resultdelong$AUC[2],resultdelong$SE[2],resultdelong$R[1,2],resultdelong$AUC[1]-resultdelong$AUC[2],resultdelong$Z,resultdelong$pvalue,result$AUC1,result$AUC2,resultboot$pvalue,nc,IC1.auc.per[1],IC1.auc.per[2],IC2.auc.per[1],IC2.auc.per[2],IC.diff.auc.per[1],IC.diff.auc.per[2],"\\\n", sep=" & ")
  #cat(resultdelong$AUC[1],resultdelong$SE[1],resultdelong$AUC[2],resultdelong$SE[2],resultdelong$R[1,2],resultdelong$AUC[1]-resultdelong$AUC[2],resultdelong$Z,resultdelong$pvalue,result$AUC1,result$AUC2,resultboot$pvalue,nc,IC1.auc.per[1],IC1.auc.per[2],IC2.auc.per[1],IC2.auc.per[2],IC.diff.auc.per[1],IC.diff.auc.per[2],"\\\n", sep=" & ",file=paste(nomeg,"Results_zang.txt"),append=TRUE)
  
  resultlist=list(Area1=resultdelong$AUC[1],SE1=resultdelong$SE[1],Area2=resultdelong$AUC[2],SE2=resultdelong$SE[2],CorrCoef=resultdelong$R[1,2],diff=(resultdelong$AUC[1]-resultdelong$AUC[2]),zstats=resultdelong$Z,pvalue1=resultdelong$pvalue,TrapArea1=result$AUC1,TrapArea2=result$AUC2,bootpvalue=resultboot$pvalue,nCross=nc,ICLB1=IC1.auc.per[1],ICUB1=IC1.auc.per[2],ICLB2=IC2.auc.per[1],ICUB2=IC2.auc.per[2],ICLBDiff=IC.diff.auc.per[1],ICUBDiff=IC.diff.auc.per[2])
  return(resultlist)
  }
