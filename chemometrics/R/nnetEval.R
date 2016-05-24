nnetEval <-
function(X,grp,train,kfold=10,decay=seq(0,10,by=1),
    size=30,maxit=100,plotit=TRUE,legend=TRUE,legpos="bottomright",...){
#
# EVALUATION for Artificial Neural Networks (ANN) by cross-validation
# change either weight decay or number of hidden units, but not both!
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

mkTable <- function(pred,tab,grplev){
  predf=factor(pred,labels=grplev[sort(unique(pred))])
  tabf=matrix(0,ncol=length(grplev),nrow=length(grplev))
  dimnames(tabf)=list(grplev,grplev)
  tabf[,levels(predf)] <- tab
  miscl <- 1-sum(diag(tabf))/sum(tabf)
  list(miscl=miscl,tab=tabf)
}

# main routine
  ntrain=length(train)
  if (min(length(decay),length(size))>1){
      stop("Either decay or size can be vectors, but not both!")
  }
  if (max(length(decay),length(size))==1){
    stop("Either decay or size must be vectors!")
  }
  else {
    lvary=max(length(decay),length(size))
    varying <- c("decay","size")[which.max(c(length(decay),length(size)))]
  }
  trainerr=rep(NA,lvary)
  testerr=rep(NA,lvary)
  cvMean=rep(NA,lvary)
  cvSe=rep(NA,lvary)
  cverr=matrix(NA,nrow=kfold,ncol=lvary)
#  require(nnet)
  for (j in 1:lvary){
    if (varying=="size"){
      resnnet=nnet(X[train,],class.ind(grp[train]),size=size[j],
                   decay=decay,maxit=maxit,trace=FALSE)
    }
    else {
      resnnet=nnet(X[train,],class.ind(grp[train]),size=size,
                   decay=decay[j],maxit=maxit,trace=FALSE)
    }
    respred=predict(resnnet,X[-train,])
    pred=apply(respred,1,which.max)
    tab=table(grp[-train],pred)
    testerr[j] <- mkTable(pred,tab,levels(grp))$miscl # test error
    respred=predict(resnnet,X[train,])
    pred=apply(respred,1,which.max)
    tab=table(grp[train],pred)
    trainerr[j] <- mkTable(pred,tab,levels(grp))$miscl # training error

    splt <- rep(1:kfold,length=ntrain)
    spltr <- sample(splt,ntrain)
    pred <- rep(NA,ntrain)
    for (i in 1:kfold){
      if (varying=="size"){
        resnnet=nnet(X[train[spltr!=i],],class.ind(grp[train[spltr!=i]]),
                   size=size[j],decay=decay,maxit=maxit,trace=FALSE)
      }
      else {
        resnnet=nnet(X[train[spltr!=i],],class.ind(grp[train[spltr!=i]]),
                   size=size,decay=decay[j],maxit=maxit,trace=FALSE)
      }
      respred=predict(resnnet,X[train[spltr==i],])
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
    if (varying=="size"){
      labx="Number of hidden units"
      labelsx=size
      plottitle=paste("Weight decay =",decay)
    }
    else {
      labx="Weight decay"
      labelsx=decay
      plottitle=paste(size,"hidden units")
    }
      
    plot(vvec,trainerr,ylim=c(0,ymax),xlab=labx,
       ylab="Missclassification error",cex.lab=1.2,type="l",lty=2,xaxt="n",...)
    title(plottitle)
    axis(1,at=vvec,labels=labelsx)
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
       cverr=cverr,decay=decay,size=size)
}

