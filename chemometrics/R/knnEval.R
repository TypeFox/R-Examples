knnEval <-
function(X,grp,train,kfold=10,knnvec=seq(2,20,by=2),
    plotit=TRUE,legend=TRUE,legpos="bottomright",...){
#
# EVALUATION for k-Nearest-Neighbors (kNN) by cross-validation
#
# subroutine for misclassification rates
evalSEfac <- function(pred,grptrain,spltr,grplev){
  kfold=max(spltr)
  k=length(grplev)
  misscli=rep(NA,k)
  for (i in 1:kfold){
    tab=table(grptrain[spltr==i],pred[spltr==i])
    misscli[i]=1-sum(diag(tab))/sum(tab)
  }
  list(mean=mean(misscli),se=sd(misscli)/sqrt(kfold),all=misscli)
}

#require(class)

# main routine
  ntrain=length(train)
  lknnvec=length(knnvec)
  trainerr=rep(NA,lknnvec)
  testerr=rep(NA,lknnvec)
  cvMean=rep(NA,lknnvec)
  cvSe=rep(NA,lknnvec)
  cverr=matrix(NA,nrow=kfold,ncol=lknnvec)
#  require(class)
  for (j in 1:lknnvec){
    pred=knn(X[train,],X[-train,],factor(grp[train]),k=knnvec[j])
    tab=table(grp[-train],pred)
    testerr[j] <- 1-sum(diag(tab))/sum(tab) # test error
    pred=knn(X[train,],X[train,],factor(grp[train]),k=knnvec[j])
    tab=table(grp[train],pred)
    trainerr[j] <- 1-sum(diag(tab))/sum(tab) # training error

    splt <- rep(1:kfold,length=ntrain)
    spltr <- sample(splt,ntrain)
    pred <- factor(rep(NA,ntrain),levels=levels(grp))
    for (i in 1:kfold){
      pred[spltr==i]=knn(X[train[spltr!=i],],X[train[spltr==i],],
               factor(grp[train[spltr!=i]]),k=knnvec[j])
    }
    resi=evalSEfac(pred,grp[train],spltr,levels(grp))
    cverr[,j] <- resi$all
    cvMean[j] <- resi$mean
    cvSe[j] <- resi$se
  }

  if (plotit){
    ymax=max(trainerr,testerr,cvMean+cvSe)
    vknnvec=seq(1,lknnvec)
    plot(vknnvec,trainerr,ylim=c(0,ymax),xlab="Number of nearest neighbors",
       ylab="Missclassification error",cex.lab=1.2,type="l",lty=2,xaxt="n",...)
    axis(1,at=vknnvec,labels=knnvec)
    points(vknnvec,trainerr,pch=4)
    lines(vknnvec,testerr,lty=1,lwd=1.3)
    points(vknnvec,testerr,pch=1)
    lines(vknnvec,cvMean,lty=1)
    points(vknnvec,cvMean,pch=16)
    for (i in 1:lknnvec){
      segments(vknnvec[i],cvMean[i]-cvSe[i],vknnvec[i],cvMean[i]+cvSe[i])
      segments(vknnvec[i]-0.2,cvMean[i]-cvSe[i],vknnvec[i]+0.2,cvMean[i]-cvSe[i])
      segments(vknnvec[i]-0.2,cvMean[i]+cvSe[i],vknnvec[i]+0.2,cvMean[i]+cvSe[i])
    }
    abline(h=min(cvMean)+cvSe[which.min(cvMean)],lty=3,lwd=1.2)      
    if (legend){
      legend(legpos,c("Test error","CV error","Training error"),
      lty=c(1,1,2),lwd=c(1.3,1,1),pch=c(1,16,4))
    }
  }
  list(trainerr=trainerr,testerr=testerr,cvMean=cvMean,cvSe=cvSe,
       cverr=cverr,knnvec=knnvec)
}

