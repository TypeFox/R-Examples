treeEval <-
function(X,grp,train,kfold=10,cp=seq(0.01,0.1,by=0.01),
    plotit=TRUE,legend=TRUE,legpos="bottomright",...){
#
# EVALUATION for Classification trees by cross-validation
#
# subroutine for misclassification rates
evalSE <- function(pred,grptrain,spltr,grplev){
  kfold=max(spltr)
  k=length(grplev)
  misscli=rep(NA,k)
  for (i in 1:kfold){
    tab=table(grptrain[spltr==i],pred[spltr==i])
    misscli[i]=mkTable(pred[spltr==i],tab,grplev)$miscl
  }
  list(mean=mean(misscli),se=sd(misscli)/sqrt(kfold),all=misscli)
}

#require(rpart)

mkTable <- function(pred,tab,grplev){
  predf=factor(pred,labels=grplev[sort(unique(pred))])
  tabf=matrix(0,ncol=length(grplev),nrow=length(grplev))
  dimnames(tabf)=list(grplev,grplev)
  tabf[,levels(predf)] <- tab
  miscl <- 1-sum(diag(tabf))/sum(tabf)
  list(miscl=miscl,tab=tabf)
}

# main routine
  dat=data.frame(grp,X)
  ntrain=length(train)
  lvary=length(cp)
  trainerr=rep(NA,lvary)
  testerr=rep(NA,lvary)
  cvMean=rep(NA,lvary)
  cvSe=rep(NA,lvary)
  cverr=matrix(NA,nrow=kfold,ncol=lvary)
  for (j in 1:lvary){
    restree=rpart(grp~.,data=dat[train,],method="class",
                   control=rpart.control(cp=.00001))
    tree2=prune(restree,cp=cp[j])
    respred=predict(tree2,newdata=dat[-train,])
    pred=apply(respred,1,which.max)
    tab=table(grp[-train],pred)
    testerr[j] <- mkTable(pred,tab,levels(grp))$miscl # test error


    respred=predict(tree2,newdata=dat[train,])
    pred=apply(respred,1,which.max)
    tab=table(grp[train],pred)
    trainerr[j] <- mkTable(pred,tab,levels(grp))$miscl # training error

    splt <- rep(1:kfold,length=ntrain)
    spltr <- sample(splt,ntrain)
    pred <- rep(NA,ntrain)
    for (i in 1:kfold){
      restree=rpart(grp~.,data=dat[train[spltr!=i],],method="class",
               control=rpart.control(cp=.00001))
      tree2=prune(restree,cp=cp[j])
      respred=predict(tree2,newdata=dat[train[spltr==i],])
      pred[spltr==i]=apply(respred,1,which.max)
    }
    resi=evalSE(pred,grp[train],spltr,levels(grp))
    cverr[,j] <- resi$all
    cvMean[j] <- resi$mean
    cvSe[j] <- resi$se
  }

  if (plotit){
    ymax=max(trainerr,testerr,cvMean+cvSe)
    vvec=seq(1,lvary)
    plot(vvec,trainerr,ylim=c(0,ymax),xlab="Tree complexity parameter",
       ylab="Missclassification error",cex.lab=1.2,type="l",lty=2,xaxt="n",...)
    axis(1,at=vvec,labels=cp)
    points(vvec,trainerr,pch=4)
    lines(vvec,testerr,lty=1,lwd=1.3)
    points(vvec,testerr,pch=1)
    lines(vvec,cvMean,lty=1)
    points(vvec,cvMean,pch=16)
    for (i in 1:lvary){
      segments(vvec[i],cvMean[i]-cvSe[i],vvec[i],cvMean[i]+cvSe[i])
      segments(vvec[i]-0.2,cvMean[i]-cvSe[i],vvec[i]+0.2,cvMean[i]-cvSe[i])
      segments(vvec[i]-0.2,cvMean[i]+cvSe[i],vvec[i]+0.2,cvMean[i]+cvSe[i])
    }
    abline(h=min(cvMean)+cvSe[which.min(cvMean)],lty=3,lwd=1.2)      
    if (legend){
      legend(legpos,c("Test error","CV error","Training error"),
      lty=c(1,1,2),lwd=c(1.3,1,1),pch=c(1,16,4))
    }
  }
  list(trainerr=trainerr,testerr=testerr,cvMean=cvMean,cvSe=cvSe,
       cverr=cverr,cp=cp)
}

