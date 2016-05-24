#######=============================#####
#-- Functions for Quantile Regression --#
#######=============================#####

cosso.qr=function(Gramat,y,tau,wt,parallel,cpus)
    {
     n<- length(y)
     p<- length(wt)

     if(p<=15)  tempM <- c(seq(0.2,floor(p/3),0.5),seq(ceiling(p/3)+1,p*0.8,0.8))
     else       tempM <- c(seq(0.5,9,.75),seq(10,ceiling(p*.7),1))

     if(parallel)  cossoObj=cosso.qr.Parallel(  Gramat,y,tau,wt,2^seq(-13,-21,-1),tempM,cpus)
     else          cossoObj=cosso.qr.Sequential(Gramat,y,tau,wt,2^seq(-13,-21,-1),tempM)

     class(cossoObj)="cosso"
     return(cossoObj)
    }

cosso.qr.Parallel=function(Gramat,y,tau,wt,tempLambda,tempM,cpus)
    {
    n<- length(y)
    p<- length(wt)
    clusterObj<-makeCluster(cpus)
    clusterEvalQ(clusterObj,library(cosso))
    clusterEvalQ(clusterObj,library(quadprog))
    clusterEvalQ(clusterObj,library(Rglpk))

   ######################################################
   CVmat  <- matrix(NA,ncol=length(tempLambda),nrow=5)
   splitID<- cvsplitID(n,5)

   for(f in 1:5)
      {
       testID <- splitID[!is.na(splitID[,f]),f]
       trainID<- (1:n)[-testID]
       trainRtheta<- wsGram(Gramat[trainID,trainID,],rep(1,p)/wt^2)
       testRtheta <- wsGram(Gramat[testID ,trainID,],rep(1,p)/wt^2)

       for(l in 1:length(tempLambda))
          {
           coefhat<- kqr(y[trainID],tau,tempLambda[l],trainRtheta)
           yhat   <- testRtheta%*%coefhat$coefs+coefhat$intercept
           CVmat[f,l]<- sum(rho(tau, y[testID]-yhat))
           }
       }
   optLam<- tempLambda[which.min(apply(CVmat,2,sum))]/ifelse(length(unique(wt))==1,4,1)
   #######################################################

     L2normMat <- matrix(NA,ncol=p,nrow=length(tempM))
     tempcoefs <- clusterApplyLB(clusterObj,tempM,twostep.qr,tau=tau,y=y,Gramat=Gramat,lam0=optLam,wt=wt)
     for(m in 1:length(tempM))
        {
        coefhat=tempcoefs[[m]]
        for(j in 1:p)  L2normMat[m,j]=sqrt(mean((coefhat$theta[j]/wt[j]^2*Gramat[,,j]%*%coefhat$coefs)^2))
        }

     if(sum(L2normMat[length(tempM),]==0)>0)  # Some component remain unselected
       {
       extMgrid=seq(max(tempM)+1,p*0.85,l=5)
       extL2normMat=matrix(NA,ncol=p,nrow=length(extMgrid))
       tempcoefs=clusterApplyLB(clusterObj,extMgrid,twostep.qr,tau=tau,y=y,Gramat=Gramat,lam0=optLam,wt=wt)
       for(m in 1:length(extMgrid))
          {
          coefhat=tempcoefs[[m]]
          for(j in 1:p)  extL2normMat[m,j]=sqrt(mean((coefhat$theta[j]/wt[j]^2*Gramat[,,j]%*%coefhat$coefs)^2))
          }
       tempM=c(tempM,extMgrid)
       L2normMat=rbind(L2normMat,extL2normMat)
       }
    stopCluster(clusterObj)
    
    cossoqrobj<- list(family="Quantile",tau=tau,cpus=cpus,basis.id=1:n,tune=list(OptLam=optLam,Mgrid=c(0,tempM),L2norm=rbind(rep(0,p),L2normMat)) )
    return(cossoqrobj)
   }


cosso.qr.Sequential=function(Gramat,y,tau,wt,tempLambda,tempM)
   {
    n<- length(y)
    p<- length(wt)
   ######################################################
   CVmat  <- matrix(NA,ncol=length(tempLambda),nrow=5)
   splitID<- cvsplitID(n,5)

   for(f in 1:5)
      {
       testID <- splitID[!is.na(splitID[,f]),f]
       trainID<- (1:n)[-testID]
       trainRtheta<- wsGram(Gramat[trainID,trainID,],rep(1,p)/wt^2)
       testRtheta <- wsGram(Gramat[testID ,trainID,],rep(1,p)/wt^2)

       for(l in 1:length(tempLambda))
          {
           coefhat<- kqr(y[trainID],tau,tempLambda[l],trainRtheta)
           yhat   <- testRtheta%*%coefhat$coefs+coefhat$intercept
           CVmat[f,l]<- sum(rho(tau, y[testID]-yhat))
           }
       }
   optLam<- tempLambda[which.min(apply(CVmat,2,sum))]/ifelse(length(unique(wt))==1,4,1)
   ########################################################
     L2normMat <- matrix(NA,ncol=p,nrow=length(tempM))
     for(m in 1:length(tempM))
        {
        coefhat <- twostep.qr(tau,y,Gramat,optLam,tempM[m],wt)
        for(j in 1:p)  L2normMat[m,j]=sqrt(mean((coefhat$theta[j]/wt[j]^2*Gramat[,,j]%*%coefhat$coefs)^2))
        }

     if(sum(L2normMat[length(tempM),]==0)>0)  # Some component remain unselected
       {
       extMgrid=seq(max(tempM)+1,p*0.85,l=5)
       extL2normMat=matrix(NA,ncol=p,nrow=length(extMgrid))
       for(m in 1:length(extMgrid))
          {
          coefhat <- twostep.qr(tau,y,Gramat,optLam,extMgrid[m],wt)
          for(j in 1:p)  extL2normMat[m,j]=sqrt(mean((coefhat$theta[j]/wt[j]^2*Gramat[,,j]%*%coefhat$coefs)^2))
          }
       tempM=c(tempM,extMgrid)
       L2normMat=rbind(L2normMat,extL2normMat)
       }
    
    cossoqrobj<- list(family="Quantile",tau=tau,cpus=1,basis.id=1:n,tune=list(OptLam=optLam,Mgrid=c(0,tempM),L2norm=rbind(rep(0,p),L2normMat)) )
    return(cossoqrobj)
   }


tune.cosso.qr <- function(object,folds=5,plot.it=TRUE)
    {
     n <- length(object$y)
     p <- length(object$wt)
     #--  Tuning Grids --#
     origMgrid<- object$tune$Mgrid[-1]
     uniqueSize=unique(apply(object$tune$L2norm[-1,]>0,1,sum))
     newGrid=origMgrid[apply(object$tune$L2norm[-1,]>0,1,sum)<=ceiling(p/3)]
     for(k in 1:length(uniqueSize))   
        {
        if(uniqueSize[k]>ceiling(p/3))  newGrid=c(newGrid,  origMgrid[max(which(apply(object$tune$L2norm[-1,]>0,1,sum)==uniqueSize[k]))]  )
        }

     newGrid=c(0,sort(newGrid))
     uniqueSize=c(0,sort(uniqueSize))
     refinePt=which(uniqueSize[-1]-uniqueSize[-length(uniqueSize)]>1)
     if(length(refinePt)>0)
         {
         refinePt1=refinePt[refinePt<10]
         refinePt2=refinePt[refinePt>=10]
         extMgrid<-             as.numeric(apply(cbind(origMgrid[refinePt1],origMgrid[refinePt1+1]),1,quantile,c(.3,.6)))
         extMgrid<-  c(extMgrid,as.numeric(apply(cbind(origMgrid[refinePt2],origMgrid[refinePt2+1]),1,mean             )) )
         }
      else
         {   extMgrid<- NULL }
     cand.M <- sort(c(newGrid[-1],extMgrid))

     IDmat<- cvsplitID(n,folds)
     cvRaw  <- matrix(NA,ncol=length(cand.M),nrow=n)
     if(object$cpus>1)
        {
         clusterObj<-makeCluster(object$cpus)
         clusterEvalQ(clusterObj,library(cosso))
         clusterEvalQ(clusterObj,library(quadprog))
         clusterEvalQ(clusterObj,library(Rglpk))
        
         for(f in 1:folds)
            {
              testID <- IDmat[!is.na(IDmat[,f]),f]
              trainID<- (1:n)[-testID]
              trainGramat<- object$Kmat[trainID,trainID,]
              testGramat <- object$Kmat[testID ,trainID,]
              #--- Parallel ---#
              coefhat<- clusterApplyLB(clusterObj,cand.M,twostep.qr,tau=object$tau,y=object$y[trainID],Gramat=trainGramat,lam0=object$tune$OptLam,wt=object$wt)
              for(m in 1:length(cand.M))
                 {
                  tempObj=coefhat[[m]]
                  yhat<- tempObj$intercept+wsGram(testGramat,tempObj$theta/object$wt^2)%*%tempObj$coefs
                  cvRaw[testID,m]=rho(object$tau,object$y[testID]-yhat)
                 }
             }
        stopCluster(clusterObj)   
        }
     else   
        {
        for(f in 1:folds)
            {
              testID <- IDmat[!is.na(IDmat[,f]),f]
              trainID<- (1:n)[-testID]
              trainGramat<- object$Kmat[trainID,trainID,]
              testGramat <- object$Kmat[testID ,trainID,]

              for(m in 1:length(cand.M))
                 {
                  tempObj=twostep.qr(object$tau,object$y[trainID],trainGramat,object$tune$OptLam,cand.M[m],object$wt)
                  yhat<- tempObj$intercept+wsGram(testGramat,tempObj$theta/object$wt^2)%*%tempObj$coefs
                  cvRaw[testID,m]=rho(object$tau,object$y[testID]-yhat)
                 }
             }
        }
     cvm =apply(cvRaw,2,mean)
     cvsd=sqrt(apply(scale(cvRaw, cvm, FALSE)^2, 2,mean)/n)

     locMinid <- which((cvm[-c(length(cvm),length(cvm)-1)]> cvm[-c(1,length(cvm))])*(cvm[-c(1,length(cvm))]<cvm[-c(1:2)])==TRUE)+1
     locMaxid <- which((cvm[-c(length(cvm),length(cvm)-1)]< cvm[-c(1,length(cvm))])*(cvm[-c(1,length(cvm))]>cvm[-c(1:2)])==TRUE)+1
     locMinid <- locMinid[locMinid<ifelse(length(locMaxid)>0,max(locMaxid),length(cvm))]
     opt.M=cand.M[which.min(cvm)]
     if(length(locMinid)>0)   opt.M=cand.M[ locMinid[which.min(cvm[locMinid[1:length(locMinid)]])] ]
     
     if(plot.it)
       {
       par(mfcol=c(1,2))
       plot(cand.M,cvm,ylim=c(min(cvm-cvsd),max(cvm+cvsd)),type="b",col=2,pch=16,lwd=1.5,xlab="M",ylab="Cross-Validated Check Error")
       for(m in 1:length(cand.M))  segments(cand.M[m],cvm[m]-cvsd[m],cand.M[m],cvm[m]+cvsd[m],col=grey(0.6))
       abline(v=opt.M,lty=2,col=2);axis(3,opt.M)
       matplot(object$tune$Mgrid,object$tune$L2norm,type="l",lty=1,col=c(1,rainbow(length(object$wt)-1)),xlab="M",ylab=expression(L[2]-norm))
       abline(v=opt.M,lty=2,col=2);axis(3,opt.M)
       axis(4,at=object$tune$L2norm[nrow(object$tune$L2norm),],labels=1:length(object$wt),cex=.3,las=2)
       }
    return(list(OptM=opt.M,M=cand.M,cvm=cvm,cvsd=cvsd))
    }

twostep.qr=function(tau,y,Gramat,lam0,mm,wt)
  {
  n <- length(y)
  p <- length(wt)
  G <- matrix(0,n,p)

  newtheta<- rep(1,p)
  #---- Step 1.2 -----#
  Rtheta  <- wsGram(Gramat,newtheta/wt^2)
  cb   <- kqr(y,tau,lam0,Rtheta,insure=TRUE)
  newc <- cb$coefs
  newb <- as.numeric(cb$intercept)
  #---- Step 2.1 -----#
  for (i in 1:p)  G[,i]=(wt[i]^(-2))*Gramat[,,i]%*%as.vector(newc)
  newtheta <- garrote.qr(x=G,y=y-newb,tau=tau,lambda0=lam0,M=mm,ct=t(newc))
  #---- Step 2.2 -----#
  Rtheta <- wsGram(Gramat,newtheta/wt^2)
  cb   <- kqr(y,tau,lam0,Rtheta,insure=TRUE)
  newc <- cb$coefs
  newb <- as.numeric(cb$intercept)

  output<-list(coefs=newc,intercept=newb,theta=newtheta,quantile=cb$quantile)
  return(output)
  }


kqr=function(y,tau,lambda,Rtheta,insure=TRUE,posconst=1e-8)
   {
    n<-length(y)
    r<-1/(2*n*lambda)
    if(insure)      Rtheta<-Rtheta+posconst*diag(1,n)
    Amat<-rbind(diag(rep(1,n)),diag(rep(-1,n)))
    Amat<-t(rbind(rep(1,n),Amat))
    b0<-c(0,rep(r*(tau-1),n),rep(-r*tau,n))
    coefhat <-drop(solve.QP(Dmat=Rtheta,dvec=y,Amat=Amat,bvec=b0,meq=1)$solution)
    fhat    <-drop(Rtheta%*%coefhat)
    intercept<-quantile(y-fhat,tau)
    output   <-list(coefs=coefhat,intercept=intercept,quantile=(fhat+intercept))
    return(output)
   }


garrote.qr=function(x,y,tau,lambda0,ct,M)
   {
   n <- nrow(x)
   p <- ncol(x)
   cvec<-c(lambda0*drop(ct%*%x)-(tau-0.5)*apply(x,2,mean),rep(1,n)/(2*n))
   mat.con<-cbind(diag(p),matrix(0,p,n))
   mat.con<-rbind(mat.con,c(-rep(1,p),rep(0,n)))
   mat.con<-rbind(mat.con,cbind(x,diag(n)))
   mat.con<-rbind(mat.con,cbind(-x,diag(n)))
   lp.out<-Rglpk_solve_LP(obj=cvec,mat=mat.con,dir=rep(">=",2*n+p+1),rhs=c(rep(0,p),-M,y,-y),max=FALSE)$solution
   theta<-lp.out[1:p]
   theta[theta<1e-5]<-0
   return(theta)
   }


rho<- function(tau,r)
   {
   sol <- tau*r*(r>0)-(1-tau)*r*(r<=0)
   return(sol)
   }


SSANOVAwt.qr<- function(x,y,tau,folds=5,cpus=1)
   {
    n<- length(y)
    p<- ncol(x)

    Gram   <- bigGram(x,x)
    Rtheta <- wsGram(Gram,rep(1,p))

    cand.lam0<- 2^seq(-10,-20,by=-0.75)
    cvMat <- matrix(NA,nrow=folds,ncol=length(cand.lam0))
    IDmat <- cvsplitID(n,folds)

    if(cpus>1)
      {
      clusterObj<-makeCluster(cpus)
      clusterEvalQ(clusterObj,library(cosso))
      clusterEvalQ(clusterObj,library(quadprog))

      for(f in 1:folds)
          {
          testID <- IDmat[!is.na(IDmat[,f]),f]
          trainID<- (1:n)[-testID]
          trainRtheta<- wsGram(Gram[trainID,trainID,],rep(1,p))
          testRtheta <- wsGram(Gram[testID ,trainID,],rep(1,p))
          #--- Parallel Computing ---#
          coefhat<- clusterApplyLB(clusterObj,cand.lam0,kqr,Rtheta=trainRtheta,y=y[trainID],tau=tau)
  
          for(l in 1:length(cand.lam0))
             {
             yhat   <- as.numeric(testRtheta%*%coefhat[[l]]$coefs+coefhat[[l]]$intercept)
             cvMat[f,l]<- sum(rho(tau, y[testID]-yhat))
             }
          }
       stopCluster(clusterObj)
       }
    else
       {
       for(f in 1:folds)
          {
          testID <- IDmat[!is.na(IDmat[,f]),f]
          trainID<- (1:n)[-testID]
          trainRtheta<- wsGram(Gram[trainID,trainID,],rep(1,p))
          testRtheta <- wsGram(Gram[testID ,trainID,],rep(1,p))
  
          for(l in 1:length(cand.lam0))
             {
             coefhat<- kqr(y[trainID],tau,cand.lam0[l],trainRtheta)
             yhat   <- as.numeric(testRtheta%*%coefhat$coefs+coefhat$intercept)
             cvMat[f,l]<- sum(rho(tau, y[testID]-yhat))
             }
          }
       }   
    optLam=cand.lam0[which.min(apply(cvMat,2,sum))]

    result=kqr(y=y,tau=tau,lambda=optLam,Rtheta=Rtheta,insure=TRUE)

    kqrnorms=rep(NA,p)
    for(j in 1:p)    kqrnorms[j]=ifelse(length(unique(x[,j]))>6,sqrt(mean((Gram[,,j]%*%result$coef)^2)),diff(range(Gram[,,j]%*%result$coef)))
    return(1/kqrnorms)
   }
