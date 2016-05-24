#######=============================#####
#------ Functions for Cox PH Model -----#
#######=============================#####
cosso.Cox <- function(Gramat,time,status,wt,nbasis,basis.id,parallel,cpus)
    {
     n<- length(time)
     p<- length(wt)
 
     if(missing(nbasis) & missing(basis.id))
       {
       nbasis=min( max(35, ceiling(11 * length(time)^(2/9)))  ,  sum(status==1)-5  )
       basis.id=sort(sample(which(status==1),nbasis))
       }
     else if(!missing(nbasis) & missing(basis.id))
       {
       if(nbasis>sum(status==1)) { nbasis=sum(status==1)-5; cat("Use nbasis=",sum(status==1)-5,"\n")}
       basis.id=sort(sample(which(status==1),nbasis))
       }
     Gramat1<-Gramat[        ,basis.id,]
     Gramat2<-Gramat[basis.id,basis.id,]
     RS <- RiskSet(time,status)

     bestlam <- ACV.lambda(Gramat,time,status,1/wt^2,basis.id,RS)
     if(p<=15)  tempM <- c(seq(0.2,floor(p/3),0.5),seq(ceiling(p/3)+1,p*0.8,0.75))
     else       tempM <- c(seq(0.5,8,.75),seq(9,p*.75,1))

     if(parallel)    cossoObj=cosso.Cox.Parallel(  Gramat1,Gramat2,time,status,wt,basis.id,RS,bestlam,tempM,cpus)
     else            cossoObj=cosso.Cox.Sequential(Gramat1,Gramat2,time,status,wt,basis.id,RS,bestlam,tempM)

     class( cossoObj)="cosso"
     return(cossoObj)
    }

cosso.Cox.Parallel <- function(Gramat1,Gramat2,time,status,wt,basis.id,RS,lambda0,initMgrid,cpus)
   {
   clusterObj=makeCluster(cpus)
   clusterEvalQ(clusterObj,library(cosso))
   clusterEvalQ(clusterObj,library(quadprog))
   clusterEvalQ(clusterObj,library(glmnet))

     p=length(wt)
     n=length(time)
     L2normMat <- matrix(NA,ncol=p,nrow=length(initMgrid))
     ACVscore  <- rep(NA,length(initMgrid))

     tempcoefs <- clusterApplyLB(clusterObj,initMgrid,twostep.Cox,Gramat1=Gramat1,Gramat2=Gramat2,time=time,status=status,wt=wt,RS=RS,lambda0=lambda0)

     for(m in 1:length(initMgrid))
        {
        coefhat <- tempcoefs[[m]]
        ACVscore[m] <- PartialLik(time,status,RS,coefhat$fit)+sum(status==1)/n^2*( sum(diag(coefhat$UHU))/(n-1) - sum(coefhat$UHU)/(n^2-n) )
        for(j in 1:p)  L2normMat[m,j] <- sqrt(mean((coefhat$theta[j]/wt[j]^2*Gramat1[,,j]%*%coefhat$coefs)^2))
        }

     if(sum(L2normMat[length(initMgrid),]==0)>0)  # Some component remain unselected
       {
       extMgrid=seq(max(initMgrid)+0.5,p*0.9,l=5)
       extL2normMat <- matrix(NA,ncol=p,nrow=length(extMgrid))
       extACVscore  <- rep(NA,length(extMgrid))
       tempcoefs=clusterApplyLB(clusterObj,extMgrid,twostep.Cox,Gramat1=Gramat1,Gramat2=Gramat2,time=time,status=status,wt=wt,RS=RS,lambda0=lambda0)
       for(m in 1:length(extMgrid))
          {
          coefhat=tempcoefs[[m]]
          extACVscore[m] <- PartialLik(time,status,RS,coefhat$fit)+sum(status==1)/n^2*( sum(diag(coefhat$UHU))/(n-1) - sum(coefhat$UHU)/(n^2-n) )
          for(j in 1:p)  extL2normMat[m,j] <- sqrt(mean((coefhat$theta[j]/wt[j]^2*Gramat1[,,j]%*%coefhat$coefs)^2))
          }
       initMgrid=c(initMgrid,extMgrid)
       ACVscore <- c(ACVscore,extACVscore)
       L2normMat=rbind(L2normMat,extL2normMat)
       }
   stopCluster(clusterObj)
   cossoqrobj<- list(family="Cox",basis.id=basis.id,RiskSet=RS,cpus=cpus,tune=list(OptLam=lambda0,ACV=c(NA,ACVscore),Mgrid=c(0,initMgrid),L2norm=rbind(rep(0,p),L2normMat)) )
   return(cossoqrobj)
   }


cosso.Cox.Sequential <- function(Gramat1,Gramat2,time,status,wt,basis.id,RS,lambda0,initMgrid)
     {
     p=length(wt)
     n=length(time)
     L2normMat <- matrix(NA,ncol=p,nrow=length(initMgrid))
     ACVscore  <- rep(NA,length(initMgrid))

     for(m in 1:length(initMgrid))
        {
        coefhat <- twostep.Cox(Gramat1,Gramat2,time,status,wt,RS,lambda0,initMgrid[m])
        ACVscore[m] <- PartialLik(time,status,RS,coefhat$fit)+sum(status==1)/n^2*( sum(diag(coefhat$UHU))/(n-1) - sum(coefhat$UHU)/(n^2-n) )
        for(j in 1:p)  L2normMat[m,j] <- sqrt(mean((coefhat$theta[j]/wt[j]^2*Gramat1[,,j]%*%coefhat$coefs)^2))
        }

     if(sum(L2normMat[length(initMgrid),]==0)>0)  # Some component remain unselected
       {
       extMgrid=seq(max(initMgrid)+0.5,p*0.85,l=5)
       extL2normMat <- matrix(NA,ncol=p,nrow=length(extMgrid))
       extACVscore  <- rep(NA,length(extMgrid))
       for(m in 1:length(extMgrid))
          {
          coefhat <- twostep.Cox(Gramat1,Gramat2,time,status,wt,RS,lambda0,extMgrid[m])
          extACVscore[m] <- PartialLik(time,status,RS,coefhat$fit)+sum(status==1)/n^2*( sum(diag(coefhat$UHU))/(n-1) - sum(coefhat$UHU)/(n^2-n) )
          for(j in 1:p)  extL2normMat[m,j] <- sqrt(mean((coefhat$theta[j]/wt[j]^2*Gramat1[,,j]%*%coefhat$coefs)^2))
          }
       initMgrid=c(initMgrid,extMgrid)
       ACVscore <- c(ACVscore,extACVscore)
       L2normMat=rbind(L2normMat,extL2normMat)
       }
    cossoqrobj<- list(family="Cox",basis.id=basis.id,RiskSet=RS,cpus=1,tune=list(OptLam=lambda0,ACV=c(NA,ACVscore),Mgrid=c(0,initMgrid),L2norm=rbind(rep(0,p),L2normMat)) )
    return(cossoqrobj)
    }

   
tune.cosso.Cox <- function(object,folds=5,plot.it=TRUE)
   {
    n=nrow(object$y)
    p=length(object$wt)

   #--  Tuning Grids --#
    origMgrid  <- object$tune$Mgrid[-1]
    uniqueSize <- unique(apply(object$tune$L2norm[-1,]>0,1,sum))
    newGrid=origMgrid[apply(object$tune$L2norm[-1,]>0,1,sum)<=ceiling(p/3)]
    for(k in 1:length(uniqueSize))   
       {
       if(uniqueSize[k]>ceiling(p/3))  newGrid=c(newGrid,  origMgrid[max(which(apply(object$tune$L2norm[-1,]>0,1,sum)==uniqueSize[k]))]  )
       }

    newGrid    <- c(0,sort(newGrid))
    uniqueSize <- c(0,sort(uniqueSize))
    refinePt   <- which(uniqueSize[-1]-uniqueSize[-length(uniqueSize)]>1)
    if(length(refinePt)>0)
       {
       refinePt1 <- refinePt[refinePt<10]
       refinePt2 <- refinePt[refinePt>=10]
       extMgrid <-            as.numeric(apply(cbind(origMgrid[refinePt1],origMgrid[refinePt1+1]),1,quantile,c(.3,.6)))
       extMgrid <- c(extMgrid,as.numeric(apply(cbind(origMgrid[refinePt2],origMgrid[refinePt2+1]),1,mean             )) )
       }
    else
       {   extMgrid<- NULL }

   cand.M <- sort(c(newGrid[-1],extMgrid))
   IDmat<- cvsplitID(n,folds)
   cvMat<- matrix(NA,ncol=length(cand.M),nrow=folds)
   if(object$cpus>1) 
      {
      clusterObj=makeCluster(object$cpus)
      clusterEvalQ(clusterObj,library(cosso))
      clusterEvalQ(clusterObj,library(quadprog))
      clusterEvalQ(clusterObj,library(glmnet))
      
      for(f in 1:folds)
         {
         testID=IDmat[!is.na(IDmat[,f]),f]
         trainID=(1:n)[-testID]
         trainRS=RiskSet(object$y[trainID,"time"],object$y[trainID,"status"])
         testRS =RiskSet(object$y[testID ,"time"],object$y[testID ,"status"])

         testGramat1 =object$Kmat[ testID,object$basis.id,]

         tempcoefs <- clusterApplyLB(clusterObj,cand.M,twostep.Cox,Gramat1=object$Kmat[trainID,object$basis.id,],Gramat2=object$Kmat[object$basis.id,object$basis.id,],time=object$y[trainID,"time"],status=object$y[trainID,"status"],wt=object$wt,RS=trainRS,lambda0=object$tune$OptLam)
         for(m in 1:length(cand.M))
            {
            tempObj =tempcoefs[[m]]
            riskScore=wsGram(testGramat1,tempObj$theta/object$wt^2)%*%tempObj$coef
            cvMat[f,m] <- PartialLik(object$y[testID,"time"],object$y[testID,"status"],testRS,riskScore)
            }
         }
      stopCluster(clusterObj)
      }
   else              
      {
      for(f in 1:folds)
         {
         testID=IDmat[!is.na(IDmat[,f]),f]
         trainID=(1:n)[-testID]
         trainRS=RiskSet(object$y[trainID,"time"],object$y[trainID,"status"])
         testRS =RiskSet(object$y[testID ,"time"],object$y[testID ,"status"])

         testGramat1=object$Kmat[testID,object$basis.id,]

         for(m in 1:length(cand.M))
             {
             tempObj <- twostep.Cox(object$Kmat[trainID,object$basis.id,],object$Kmat[object$basis.id,object$basis.id,],time=object$y[trainID,"time"],object$y[trainID,"status"],object$wt,trainRS,object$tune$OptLam,cand.M[m])
             riskScore<-wsGram(testGramat1,tempObj$theta/object$wt^2)%*%tempObj$coef
             cvMat[f,m] <- PartialLik(object$y[testID,"time"],object$y[testID,"status"],testRS,riskScore)
             }
         }    
      }
   cvm = apply(cvMat,2,mean)
   cvsd= sqrt(apply(scale(cvMat, cvm, FALSE)^2, 2, mean)/(folds- 1))
        
   locMinid <- which((cvm[-c(length(cvm),length(cvm)-1)]> cvm[-c(1,length(cvm))])*(cvm[-c(1,length(cvm))]<cvm[-c(1:2)])==TRUE)+1
   locMaxid <- which((cvm[-c(length(cvm),length(cvm)-1)]< cvm[-c(1,length(cvm))])*(cvm[-c(1,length(cvm))]>cvm[-c(1:2)])==TRUE)+1
   locMaxid <- locMaxid[locMaxid>(length(cand.M)*2/3)]
   locMinid <- locMinid[locMinid<ifelse(length(locMaxid)>0,max(locMaxid),length(cvm))]
   opt.M=cand.M[which.min(cvm)]
   if(length(locMinid)>0)   opt.M=ifelse( length(locMinid)>2,cand.M[ locMinid[which.min(cvm[locMinid[1:(length(locMinid)-1)]])] ],cand.M[locMinid[ which.min(cvm[locMinid[1:max(1,length(locMinid))]])] ] )

   if(plot.it)
      {
       par(mfcol=c(1,2))
       plot(cand.M,cvm,ylim=c(min(cvm-cvsd),max(cvm+cvsd)),type="b",col=2,pch=16,lwd=1.5,xlab="M",ylab="Cross-Validated minus log-Partial Likelihood")
       for(m in 1:length(cand.M))  segments(cand.M[m],cvm[m]-cvsd[m],cand.M[m],cvm[m]+cvsd[m],col=grey(0.6))
       abline(v=opt.M,lty=2,col=2);axis(3,opt.M)
       matplot(object$tune$Mgrid,object$tune$L2norm,type="l",lty=1,col=c(1,rainbow(length(object$wt)-1)),xlab="M",ylab=expression(L[2]-norm))
       abline(v=opt.M,lty=2,col=2);axis(3,opt.M)
       axis(4,at=object$tune$L2norm[nrow(object$tune$L2norm),],labels=1:length(object$wt),cex=.3,las=2)
       }
    return(list(OptM=opt.M,M=cand.M,cvm=cvm,cvsd=cvsd))
   }


    

ACV.lambda <- function(Gramat,time,status,mscale,basis.id,RS,cand.lambda0)
   {
   if(missing(cand.lambda0)) cand.lambda0=2^seq(-10,-21,-0.75)
   cand.lambda0=sort(cand.lambda0,decreasing=TRUE)

   n=length(time)
   p=length(mscale)

   tempBasisID=basis.id
   if(length(basis.id)>30)   tempBasisID=sort(sample(basis.id,30))

   Rtheta1=wsGram(Gramat[,tempBasisID,],mscale)

   Hess.FullNumer.unScale=array(NA,dim=c(length(tempBasisID),length(tempBasisID),n))
   for(i in 1:n)  Hess.FullNumer.unScale[,,i] =Rtheta1[i,]%*%t(Rtheta1[i,])

   ACV=matrix(NA,nrow=length(cand.lambda0),ncol=2)
   for(j in 1:length(cand.lambda0))
      {
      tempCox=sspline.Cox(Gramat[,tempBasisID,],Gramat[tempBasisID,tempBasisID,],time,status,mscale,cand.lambda0[j],RS,Hess.FullNumer.unScale)
      ACV[j,1]=PartialLik(time,status,RS,tempCox$fit)
      ACV[j,2]=ACV[j,1]+sum(status==1)/n^2*( sum(diag(tempCox$UHU))/(n-1) - sum(tempCox$UHU)/(n^2-n) )
      }
   acv=ACV[,2]
   locMinid <- which((acv[-c(length(acv),length(acv)-1)]> acv[-c(1,length(acv))])*(acv[-c(1,length(acv))]<acv[-c(1:2)])==TRUE)+1
   locMaxid <- which((acv[-c(length(acv),length(acv)-1)]< acv[-c(1,length(acv))])*(acv[-c(1,length(acv))]>acv[-c(1:2)])==TRUE)+1
   locMinid <- locMinid[locMinid<ifelse(length(locMaxid)>0,max(locMaxid),length(acv))]
   opt.lambda0=cand.lambda0[which.min(acv)]
   if(length(locMinid)>0)   opt.lambda0=cand.lambda0[ locMinid[which.min(acv[locMinid[1:length(locMinid)]])] ]

   return(opt.lambda0/ifelse(length(unique(mscale))==1,4,1))
   }


#---- Solve SS-ANOVA Cox and then a Quadratic Programming ----#
twostep.Cox <- function(Gramat1,Gramat2,time,status,wt,RS,lambda0,M)
  {
  n=length(time)
  p=length(wt)
  #---- Step 1.2 ----#
  ssCox =sspline.Cox(Gramat1,Gramat2,time,status,    rep(1,p)/wt^2,lambda0,RS)
  L2norm=rep(NA,p)
  for(j in 1:p)   L2norm[j]=sqrt(mean( (wt[j]^(-2)*Gramat1[,,j]%*%ssCox$coefs)^2 ))
  init.Theta=L2norm*M/sum(L2norm)
  #---- Step 2.1 ----#
   garCox=garrote.Cox(Gramat1,Gramat2,time,status,wt,lambda0,M,ssCox$coef,init.Theta,RS)
  #---- Step 2.2 ----#
   ssCox =sspline.Cox(Gramat1,Gramat2,time,status,garCox$theta/wt^2,lambda0,RS)
  obj=c(ssCox,list(theta=garCox$theta,intercept=0))

  return(obj)
  }

#---- Solve a SS-ANOVA Cox Problem ----#
sspline.Cox <- function(Gramat1,Gramat2,time,status,mscale,lambda0,RS,Hess.FullNumer.unScale)
  {
  n=length(time)
  p=length(mscale)
  Rtheta1=wsGram(Gramat1,mscale)
  Rtheta2=wsGram(Gramat2,mscale)

  EigRtheta2=eigen(Rtheta2)
  if(min(EigRtheta2$value)<0)
     {
     Rtheta2=Rtheta2+max(1e-7,1.5*abs(min(EigRtheta2$value)))*diag(nrow(Rtheta2))
     EigRtheta2=eigen(Rtheta2)
     }
  pseudoX=Rtheta1%*%EigRtheta2$vectors%*%diag(sqrt(1/EigRtheta2$values))
  ssCox.en=glmnet(pseudoX,cbind(time=time,status=status),family="cox",lambda=c(lambda0/2,lambda0),alpha=0,standardize = FALSE)
  init.C=as.numeric( EigRtheta2$vectors%*%diag(sqrt(1/EigRtheta2$values))%*%ssCox.en$beta[,1] )

  #---- One-Step Update ----#
  f.old=Rtheta1%*%init.C
  GH=gradient.Hessian.C(init.C,Gramat1,Gramat2,time,status,mscale,lambda0,RS,Hess.FullNumer.unScale)
  new.C=as.numeric(My_solve(GH$H,GH$H%*%init.C-GH$G))


  UHU=Rtheta1%*%My_solve(GH$H,t(Rtheta1))
  ssCoxObj=list(coefs=new.C,fit=Rtheta1%*%new.C,UHU=UHU)
  return(ssCoxObj)
  }

gradient.Hessian.C=function(initC,Gramat1,Gramat2,time,status,mscale,lambda0,riskset,Hess.FullNumer.unScale)
  {
  n=length(time)
  tie.size=as.numeric( table(time[status==1]) )

  Rtheta1=wsGram(Gramat1,mscale)
  Rtheta2=wsGram(Gramat2,mscale)
  if(min(eigen(Rtheta2)$value)<0)  Rtheta2=Rtheta2+1e-8*diag(nrow(Rtheta2))
  eta=Rtheta1%*%initC

  if(missing(Hess.FullNumer.unScale))
     {
     Hess.FullNumer.unScale=array(NA,dim=c(length(initC),length(initC),n))
     for(i in 1:n)  Hess.FullNumer.unScale[,,i] =Rtheta1[i,]%*%t(Rtheta1[i,])
     }

  Grad.Term1=-t(Rtheta1)%*%status/n
  Grad.Term2=matrix(NA,ncol=ncol(riskset),nrow=length(initC))
  Grad.Term3=2*lambda0*Rtheta2%*%initC

  Grad.FullNumer=t(Rtheta1)%*%diag(as.numeric(exp(eta)))
  Grad.FullDenom=Hess.FullDenom=exp(eta)

  Hess.FullNumer =Hess.FullNumer.unScale*array( rep( exp(eta),each=length(initC)^2 ), dim=c(length(initC),length(initC),n))
  Hess.Term1=Hess.Term2=array(NA,dim=c(length(initC),length(initC),ncol(riskset)))

  k=1
  tempSum.exp.eta=sum( exp(eta[ riskset[,k] ]),na.rm=TRUE )
  temp.Gradient.numer=apply(Grad.FullNumer[, riskset[,k] ],1     ,sum,na.rm=TRUE)
  temp.Hessian.numer =apply(Hess.FullNumer[,,riskset[,k] ],c(1,2),sum,na.rm=TRUE)

  Grad.Term2[,k] =tie.size[k]*temp.Gradient.numer/tempSum.exp.eta
  Hess.Term1[,,k]=temp.Hessian.numer /tempSum.exp.eta
  Hess.Term2[,,k]=1/tie.size[k]*Grad.Term2[,k]%*%t(Grad.Term2[,k])

  for(k in 2:ncol(riskset))
     {
     excludeID=riskset[,k-1][ !riskset[,k-1]%in%riskset[,k] ]

     tempSum.exp.eta=tempSum.exp.eta-sum(exp(eta[excludeID]))
     if(length(excludeID)>1)
        {
        temp.Gradient.numer=temp.Gradient.numer-apply(Grad.FullNumer[, excludeID],1     ,sum)
        temp.Hessian.numer =temp.Hessian.numer -apply(Hess.FullNumer[,,excludeID],c(1,2),sum)
        }
     else
        {
        temp.Gradient.numer=temp.Gradient.numer-      Grad.FullNumer[, excludeID]
        temp.Hessian.numer =temp.Hessian.numer -      Hess.FullNumer[,,excludeID]
        }

      Grad.Term2[,k] =tie.size[k]*temp.Gradient.numer/tempSum.exp.eta
      Hess.Term1[,,k]=temp.Hessian.numer /tempSum.exp.eta
      Hess.Term2[,,k]=1/tie.size[k]*Grad.Term2[,k]%*%t(Grad.Term2[,k])
     }
  Grad.Term2=apply(Grad.Term2,1,sum)/n

  Gradient=Grad.Term1+Grad.Term2+Grad.Term3
  Hessian =apply(Hess.Term1,c(1,2),sum)/n-apply(Hess.Term2,c(1,2),sum)/n+2*lambda0*Rtheta2

  return(list(Gradient=Gradient,Hessian=Hessian))
  }


#---- Solve a Quadratic Programming Problem ------#
garrote.Cox <- function(Gramat1,Gramat2,time,status,wt,lambda0,M,init.C,init.Theta,RS)
  {
  n=length(time)
  p=length(wt)

  if(missing(init.Theta))
     {
     L2norm=rep(NA,p)
     for(j in 1:p)   L2norm[j]=sqrt(mean((1/wt[j]^2*Gramat1[,,j]%*%init.C)^2))
     init.Theta=L2norm*M/sum(L2norm)
     }

  G1=matrix(NA,ncol=p,nrow=n)
  G2=matrix(NA,ncol=p,nrow=dim(Gramat2)[1])
  for(j in 1:p)     
    {
    G1[,j]=1/wt[j]^2*Gramat1[,,j]%*%init.C
    G2[,j]=1/wt[j]^2*Gramat2[,,j]%*%init.C
    }

  Hess.FullNumer.unScale=array(NA,dim=c(length(init.Theta),length(init.Theta),n))
  for(i in 1:n)  Hess.FullNumer.unScale[,,i] =G1[i,]%*%t(G1[i,])

  loop=0
  iter.diff=Inf
  old.Theta=init.Theta
  while(loop<15 & iter.diff>1e-4)
      {
      loop=loop+1
      GH=gradient.Hessian.Theta(old.Theta,init.C,G1,G2,lambda0,M,time,status,RS,Hess.FullNumer.unScale)
      if(min(eigen(GH$H)$value)<0)   GH$H=GH$H+max( 1e-7,1.5*abs(min(eigen(GH$H)$value)) )*diag(length(init.Theta))
      dvec=-(GH$G-GH$H%*%old.Theta)
      Amat=t(rbind(diag(p),rep(-1,p)))
      bvec=c(rep(0,p),-M)
      new.Theta=My_solve.QP(GH$H,dvec,Amat,bvec)
      new.Theta[new.Theta<1e-7]=0
      iter.diff=mean(abs(new.Theta-old.Theta))
      old.Theta=new.Theta
      }
  return(list(coefs=init.C,theta=new.Theta))
  }

gradient.Hessian.Theta=function(initTheta,initC,G1,G2,lambda0,M,time,status,riskset,Hess.FullNumer.unScale)
  {
  n=length(time)
  p=length(initTheta)
  tie.size=as.numeric( table(time[status==1]) )
  eta=G1%*%initTheta

  Grad.Term1=-t(G1)%*%status/n
  Grad.Term2=matrix(NA,ncol=ncol(riskset),nrow=p)
  Grad.Term3=lambda0*t(G2)%*%initC

  Grad.FullNumer=t(G1)%*%diag(as.numeric(exp(eta)))
  Grad.FullDenom=Hess.FullDenom=exp(eta)

  Hess.FullNumer =Hess.FullNumer.unScale*array( rep( exp(eta),each=p^2 ), dim=c(p,p,n))
  Hess.Term1=Hess.Term2=array(NA,dim=c(p,p,ncol(riskset)))

  k=1
     tempSum.exp.eta=sum( exp( eta[ riskset[,k] ] ) ,na.rm=TRUE)
     tempGradient.numer=apply(Grad.FullNumer[,  riskset[,k] ],1     ,sum,na.rm=TRUE)
     tempHessian.numer =apply(Hess.FullNumer[,, riskset[,k] ],c(1,2),sum,na.rm=TRUE)

     Grad.Term2[,k] =tie.size[k]*tempGradient.numer/tempSum.exp.eta
     Hess.Term1[,,k]=            tempHessian.numer /tempSum.exp.eta
     Hess.Term2[,,k]=1/tie.size[k]*Grad.Term2[,k]%*%t(Grad.Term2[,k])

  for(k in 2:ncol(riskset))
     {
     excludeID=riskset[,k-1][! riskset[,k-1]%in%riskset[,k] ]
     tempSum.exp.eta=tempSum.exp.eta-sum( exp(eta[excludeID]) )

     if(length(excludeID)>1)
        {
        tempGradient.numer=tempGradient.numer-apply(Grad.FullNumer[, excludeID],1     ,sum)
        tempHessian.numer =tempHessian.numer -apply(Hess.FullNumer[,,excludeID],c(1,2),sum)
        }
     else
        {
        tempGradient.numer=tempGradient.numer-Grad.FullNumer[, excludeID]
        tempHessian.numer =tempHessian.numer -Hess.FullNumer[,,excludeID]
        }
     Grad.Term2[,k] =tie.size[k]*tempGradient.numer/tempSum.exp.eta
     Hess.Term1[,,k]=            tempHessian.numer /tempSum.exp.eta
     Hess.Term2[,,k]=1/tie.size[k]*Grad.Term2[,k]%*%t(Grad.Term2[,k])
     }
  Grad.Term2=apply(Grad.Term2,1,sum)/n

  Gradient=Grad.Term1+Grad.Term2+Grad.Term3
  Hessian =apply(Hess.Term1,c(1,2),sum)/n-apply(Hess.Term2,c(1,2),sum)/n

  return(list(Gradient=Gradient,Hessian=Hessian))
  }

SSANOVAwt.Cox <- function(x,time,status,mscale=rep(1,ncol(x)),nbasis,basis.id)
  {
    n <- length(time)
    p <- ncol(x)
    Gramat <- bigGram(x,x)
    basis.id <- sort(sample(which(status==1),min(30+ceiling(n/100),sum(status==1)-5)))
  
    if( missing(nbasis) &  missing(basis.id))
      {
      nbasis=max(35, ceiling(10 * length(time)^(2/9)))
      basis.id=sort(sample(which(status==1),min(nbasis,sum(status==1)-5)))
      }
    if( missing(nbasis) & !missing(basis.id))     nbasis <- length(basis.id)
    if(!missing(nbasis) &  missing(basis.id))     basis.id <- sort( sample(1:n,min(nbasis,sum(status==1)-5)) )


    RS <- RiskSet(time,status)
    optLambda <- ACV.lambda(Gramat,time,status,mscale,basis.id,RS,2^seq(-8,-21,-1))
    coxObj <- sspline.Cox(Gramat[,basis.id,],Gramat[basis.id,basis.id,],time,status,mscale,optLambda,RS)

    L2norm=rep(NA,p)
    for(j in 1:p)    L2norm[j]=ifelse(length(unique(x[,j]))>6,sqrt(mean( (mscale[j]*Gramat[,basis.id,j]%*%coxObj$coefs)^2 )),diff(range((mscale[j]*Gramat[,basis.id,j]%*%coxObj$coefs))))

  return(1/L2norm)
  }


PartialLik=function(time,status,RS,fhat)
   {
   pl=rep(NA,ncol(RS))
   eventtime=unique(time[status==1])
   tie.size=as.numeric(table(time[status==1]))
   for(k in 1:ncol(RS))
      {
      failid=which(time==eventtime[k])
      pl[k]=sum(fhat[failid])-tie.size[k]*log( sum(exp( fhat[RS[,k]]),na.rm=T ) )
      }
   return(-sum(pl)/length(time))
   }



RiskSet <-  function(time,status)
  {
  eventTime=sort(unique(time[status==1]))
  RiskSet=matrix(NA,ncol=length(eventTime),nrow=length(time))
  for(k in 1:length(eventTime))
     {
     risk.id=which(time>=eventTime[k])     
     RiskSet[risk.id,k]=risk.id
     }
  return(RiskSet)
  }
