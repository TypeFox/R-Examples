svmEval <-
function(X,grp,train,kfold=10,gamvec=seq(0,10,by=1),
    kernel="radial",degree=3,plotit=TRUE,legend=TRUE,legpos="bottomright",...){
#
# EVALUATION for Support Vector Machines (SVM) by cross-validation
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

# main routine
  ntrain=length(train)
  lgamvec=length(gamvec)
  trainerr=rep(NA,lgamvec)
  testerr=rep(NA,lgamvec)
  cvMean=rep(NA,lgamvec)
  cvSe=rep(NA,lgamvec)
  cverr=matrix(NA,nrow=kfold,ncol=lgamvec)
#  require(e1071)
  for (j in 1:lgamvec){
    ressvm=svm(X[train,],factor(grp[train]),kernel=kernel,degree=3,gamma=gamvec[j])
    pred=predict(ressvm,X[-train,])
    tab=table(grp[-train],pred)
    testerr[j] <- 1-sum(diag(tab))/sum(tab) # test error
    pred=predict(ressvm,X[train,])
    tab=table(grp[train],pred)
    trainerr[j] <- 1-sum(diag(tab))/sum(tab) # training error

    splt <- rep(1:kfold,length=ntrain)
    spltr <- sample(splt,ntrain)
    pred <- factor(rep(NA,ntrain),levels=levels(grp))
    for (i in 1:kfold){
      res <- svm(X[train[spltr!=i],],factor(grp[train[spltr!=i]]),
                 kernel=kernel,degree=3,gamma=gamvec[j])
      pred[spltr==i] <- predict(res,X[train[spltr==i],])
    }
    resi=evalSEfac(pred,grp[train],spltr,levels(grp))
    cverr[,j] <- resi$all
    cvMean[j] <- resi$mean
    cvSe[j] <- resi$se
  }

  if (plotit){
    ymax=max(trainerr,testerr,cvMean+cvSe)
    vgamvec=seq(1,lgamvec)
    plot(vgamvec,trainerr,ylim=c(0,ymax),xlab="Gamma",
       ylab="Missclassification error",cex.lab=1.2,type="l",lty=2,xaxt="n",...)
    axis(1,at=vgamvec,labels=gamvec)
    points(vgamvec,trainerr,pch=4)
    lines(vgamvec,testerr,lty=1,lwd=1.3)
    points(vgamvec,testerr,pch=1)
    lines(vgamvec,cvMean,lty=1)
    points(vgamvec,cvMean,pch=16)
    for (i in 1:lgamvec){
      segments(vgamvec[i],cvMean[i]-cvSe[i],vgamvec[i],cvMean[i]+cvSe[i])
      segments(vgamvec[i]-0.2,cvMean[i]-cvSe[i],vgamvec[i]+0.2,cvMean[i]-cvSe[i])
      segments(vgamvec[i]-0.2,cvMean[i]+cvSe[i],vgamvec[i]+0.2,cvMean[i]+cvSe[i])
    }
    abline(h=min(cvMean)+cvSe[which.min(cvMean)],lty=3,lwd=1.2)      
    if (legend){
      legend(legpos,c("Test error","CV error","Training error"),
      lty=c(1,1,2),lwd=c(1.3,1,1),pch=c(1,16,4))
    }
  }
  list(trainerr=trainerr,testerr=testerr,cvMean=cvMean,cvSe=cvSe,
       cverr=cverr,gamvec=gamvec)
}

