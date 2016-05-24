vblpcmroc<-function(v.params, NUM=1e2)
  {
  N<-v.params$N
  Y<-v.params$Y
  delete<-seq(from=1, to=N*N, by=(N+1)) 
  probs<-predict.vblpcm(v.params)
  
  true_pos<-rep(NaN,NUM)
  #theta<-seq(1,0,length.out=NUM)
  theta<-quantile(probs,probs=seq(1,0,l=NUM))
  false_pos<-seq(0,1,length.out=NUM)
  tmpfunc<-function(x) abs(sum((probs[-delete]>x)[Y[-delete]==0],na.rm=1)-false_pos[i]*sum(Y==0,na.rm=1))
  
  
  for (i in 1:length(theta))
    {
    false_pos[i]=sum((probs[-delete]>theta[i])[Y[-delete]==0],na.rm=1)/sum(Y[-delete]==0,na.rm=1)
    #theta[i]=optimize(tmpfunc, interval=c(0,1))$minimum
    #theta[i]=optim(false_pos[i], tmpfunc, method="SANN", control=list(maxit=1e3))$par
    true_pos[i]=sum((probs[-delete]>theta[i])[Y[-delete]==1],na.rm=1)/sum(Y[-delete]==1,na.rm=1)
    }
  tl<-length(true_pos)
  AUC=sum(diff(false_pos)*apply(matrix(c(true_pos[1:(tl-1)],true_pos[2:tl]),ncol=2),1,mean))
  cat("AUC = ", AUC, "\n")
  plot(false_pos,true_pos,t='l',col=2,xlab="false positive rate", ylab="true positive rate",
       main=paste("AUC = ",round(AUC,3),sep=""))
  points(c(0,1),c(0,1),t='l')
  return(AUC)
  }
