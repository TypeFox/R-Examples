tri.basic <-
function(tri,M.exp,E.exp,T.exp,index=1,N=0.25){
  M.exp<-as.matrix(M.exp)
  E.exp<-as.matrix(E.exp)
  T.exp<-as.matrix(T.exp)
  
  m<-matrix(nrow=0,ncol=7)
  colnames(m)<-c("pTFtar","pTF","ptar","rlow","rhigh","plow","phigh")
  for(ind in index){
    lnc<-as.character(tri[ind,1])
    TF<-as.character(tri[ind,2])
    tar<-as.character(tri[ind,3])
    
    exp.threshold<-quantile(M.exp[lnc,],probs=N)
    indexlow<-M.exp[lnc,]<exp.threshold
    tarlow<-T.exp[tar,indexlow]
    TFlow<-E.exp[TF,indexlow]
    corl<-cor.test(T.exp[tar,indexlow],E.exp[TF,indexlow])
    plow<-corl$p.value
    rlow<-corl$estimate
    exp.threshold<-quantile(M.exp[lnc,],probs=1-N)
    indexhigh<-M.exp[lnc,]>exp.threshold
    tarhigh<-T.exp[tar,indexhigh]
    TFhigh<-E.exp[TF,indexhigh]
    corl<-cor.test(T.exp[tar,indexhigh],E.exp[TF,indexhigh])
    phigh<-corl$p.value
    rhigh<-corl$estimate
    
    pTFtar<-summary(lm(T.exp[tar,]~E.exp[TF,]))$coefficients[8]
    pTF<-t.test(E.exp[TF,indexlow],E.exp[TF,indexhigh])$p.value
    ptar<-t.test(T.exp[tar,indexlow],T.exp[tar,indexhigh])$p.value
    tmpm<-c(pTFtar,pTF,ptar,rlow,rhigh,plow,phigh)
    m<-rbind(m,tmpm)
    pdf(file=paste("basic_",paste(ind,lnc,TF,tar,sep="_"),".pdf",sep=""))
    
    laym<-matrix(2,6,6)
    laym[5:6,1:2]<-0
    laym[1:4,1:2]<-1
    laym[5:6,3:6]<-3
    layout(laym)
    options(digits=3)

    par(mar=c(0,4.1,4.1,0))
    boxplot(tarlow,tarhigh,main="",col=c("green4","red2"),names=c("LOW","HIGH"),xlab=round(ptar,3),boxwex=c(0.5,0.5),ylim=c(min(tarlow,tarhigh),max(tarlow,tarhigh)),xaxt='n')
    axis(side=3,at=c(1,2),labels=c("LOW","HIGH"),font=16)
    par(mar=c(0,0,4.1,4.1))
    plot(tarlow~TFlow,main="",xlab="",ylab="",ylim=c(min(tarlow,tarhigh),max(tarlow,tarhigh)),xlim=c(min(TFlow,TFhigh),max(TFlow,TFhigh)),xaxt='n',yaxt='n',pch=21,bg='green4')
    low_model=lm(tarlow~TFlow)
    abline(low_model,lty=1,col='green4')
    points(tarhigh~TFhigh,pch=21,bg='red2')
    high_model=lm(tarhigh~TFhigh)
    abline(high_model,lty=1,col='red2')
    par(mar=c(5.1,0,0,4.1))
    boxplot(TFlow,TFhigh,main="",col=c("green4","red2"),names=c("LOW","HIGH"),mai=0,boxwex=c(0.5,0.5),ylim=c(min(TFlow,TFhigh),max(TFlow,TFhigh)),horizontal=T,yaxt='n')
    axis(side=4,at=c(1,2),labels=c("LOW","HIGH"),font=16)
    
    dev.off()
  }#for ind(the index of triplet) 
  rownames(m)<-index
  colnames(m)<-c("P_effector_target","P_effector","P_target","R_low","R_high","P_low","P_high")
  m<-as.data.frame(m,stringsAsFactors=F)
  cat("The plot(s) is on the working directory!\n\n")
  return(m)
}
