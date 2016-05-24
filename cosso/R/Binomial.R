cosso.Binomial <- function(Gramat,y,wt,nbasis,basis.id)
  {
  n <- length(y)
  p <- length(wt)
  if( missing(nbasis) &  missing(basis.id))
    {
    nbasis=max(40, ceiling(12*n^(2/9)))
    basis.id=sort(sample(1:n,nbasis))
    }
  if(!missing(nbasis) &  missing(basis.id))     basis.id <- sort(sample(1:n,nbasis))

  Gramat1 <- Gramat[        ,basis.id,]
  Gramat2 <- Gramat[basis.id,basis.id,]

  optLambda0=cvlam.logistic(Gramat1,Gramat2,y,1/wt^2,5,2^(seq(-12,-22,-1)),TRUE)
  L2normMat <- Mgrid <- NULL
  tempM=0.2
  tempTheta=tempL2norm=rep(0,p)
  loop=0
  while(sum(tempTheta>1e-7)<p  &  loop<=ifelse(p<=15,floor(2*p),p) )
       {
       loop=loop+1
       Mgrid <- c(Mgrid,tempM)
       tempcoefhat <- twostep.Binomial(Gramat1,Gramat2,y,wt,optLambda0,tempM)
       for(j in 1:p)  tempL2norm[j] <- ifelse( length(unique(as.numeric(Gramat[,,j])))>6 , sqrt(mean((tempcoefhat$theta[j]/wt[j]^2*Gramat1[,,j]%*%tempcoefhat$coefs)^2)) , diff(range(tempcoefhat$theta[j]/wt[j]^2*Gramat1[,,j]%*%tempcoefhat$coefs)) )
       tempTheta <- tempcoefhat$theta
       L2normMat <- c(L2normMat,tempL2norm)

       if(loop<10)                   tempM <- tempM+0.25
       else if(10<=loop & loop<16)   tempM <- tempM+0.5
       else if(16<=loop & loop<20)   tempM <- tempM+1
       else                          tempM <- tempM+2
       }
  L2normMat=matrix(L2normMat,ncol=p,byrow=TRUE)
  L2normMat=rbind(rep(0,p),L2normMat)
  Mgrid=c(0,Mgrid)

  cossoObj<- list(family="Binomial",basis.id=basis.id,tune=list(OptLam=optLambda0,Mgrid=Mgrid,L2norm=L2normMat) )
  return(cossoObj)
  }

tune.cosso.Binomial<- function(object,folds,plot.it=TRUE)
   {
   #--  Tuning Grids --#
   origMgrid <- object$tune$Mgrid[-1]
   uniqueSize<- unique(apply(object$tune$L2norm[-1,]>0,1,sum))
   newGrid   <- origMgrid[apply(object$tune$L2norm[-1,]>0,1,sum)<=floor(length(object$wt)/2)]
   for(k in 1:length(uniqueSize))   
      {
      if(uniqueSize[k]>floor(length(object$wt)/2))  newGrid=c(newGrid,  origMgrid[max(which(apply(object$tune$L2norm[-1,]>0,1,sum)==uniqueSize[k]))]  )
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
   IDmat  <- cvsplitID(length(object$y),folds)
   cvRaw  <- matrix(NA,ncol=length(cand.M),nrow=length(object$y))
   for(f in 1:folds)
      {
      testID=IDmat[!is.na(IDmat[,f]),f]
      trainID=(1:length(object$y))[-testID]
      trainGramat1=object$Kmat[trainID,object$basis.id,]
      testGramat1 =object$Kmat[ testID,object$basis.id,]
      for(m in 1:length(cand.M))
         {
         tempObj =twostep.Binomial(trainGramat1,object$Kmat[object$basis.id,object$basis.id,],object$y[trainID],object$wt,object$tune$OptLam,cand.M[m])
         test.fit=tempObj$intercept+wsGram(testGramat1,tempObj$theta/object$wt^2)%*%tempObj$coefs
         cvRaw[testID,m]=log(1+exp(test.fit))-object$y[testID]*test.fit
         }
      }
    cvm <- apply(cvRaw,2,mean)
    cvsd<- sqrt(apply(scale(cvRaw, cvm, FALSE)^2, 2,mean)/length(object$y))

  locMinid <- which((cvm[-c(length(cvm),length(cvm)-1)]> cvm[-c(1,length(cvm))])*(cvm[-c(1,length(cvm))]<cvm[-c(1:2)])==TRUE)+1
  locMaxid <- which((cvm[-c(length(cvm),length(cvm)-1)]< cvm[-c(1,length(cvm))])*(cvm[-c(1,length(cvm))]>cvm[-c(1:2)])==TRUE)+1
  locMaxid <- locMaxid[locMaxid>(length(cand.M)*2/3)]
  locMinid <- locMinid[locMinid<ifelse(length(locMaxid)>0,max(locMaxid),length(cvm))]
  opt.M=cand.M[which.min(cvm)]
  if(length(locMinid)>0)   opt.M=ifelse( length(locMinid)>2,cand.M[ locMinid[which.min(cvm[locMinid[1:(length(locMinid)-1)]])] ],cand.M[locMinid[ which.min(cvm[locMinid[1:max(1,length(locMinid))]])] ] )

  if(plot.it)
     {
      par(mfcol=c(1,2))
      plot(cand.M,cvm,ylim=c(min(cvm-cvsd),max(cvm+cvsd)),type="b",col=2,pch=16,lwd=1.5,xlab="M",ylab="Cross-Validated minus log-likelihood")
      for(m in 1:length(cand.M))  segments(cand.M[m],cvm[m]-cvsd[m],cand.M[m],cvm[m]+cvsd[m],col=grey(0.6))
      abline(v=opt.M,lty=2,col=2);axis(3,opt.M)
      matplot(object$tune$Mgrid,object$tune$L2norm,type="l",lty=1,col=c(1,rainbow(length(object$wt)-1)),xlab="M",ylab=expression(L[2]-norm))
      abline(v=opt.M,lty=2,col=2);axis(3,opt.M)
      axis(4,at=object$tune$L2norm[nrow(object$tune$L2norm),],labels=1:length(object$wt),cex=.3,las=2)
     }
    return(list(OptM=opt.M,M=cand.M,cvm=cvm,cvsd=cvsd))
  }



cvlam.logistic<-function(Gramat1,Gramat2,y,mscale,folds,cand.lambda,restrict=FALSE)
  {
  cand.lambda=sort(cand.lambda,decreasing=T)

  Rtheta1=wsGram(Gramat1,mscale)
  Rtheta2=wsGram(Gramat2,mscale)
  EigRtheta2=eigen(Rtheta2);
  loop=0
  while(min(EigRtheta2$values)<0 & loop<10)
     {
     loop=loop+1
     Rtheta2=Rtheta2+1e-8*diag(nrow(Rtheta2))
     EigRtheta2=eigen(Rtheta2);
     }
   if(loop==10)   EigRtheta2$values[EigRtheta2$values<0]=1e-8
   pseudoX=Rtheta1%*%EigRtheta2$vectors%*%diag(sqrt(1/EigRtheta2$values))

   cvObj=cv.glmnet(x=pseudoX, y=y,lambda=cand.lambda,nfolds=5,alpha=0,standardize = FALSE,family="bin")

   cvs=cvObj$cvm
   opt.lambda=cand.lambda[which.min(cvs)]
   locMinid <- which((cvs[-c(length(cvs),length(cvs)-1)]> cvs[-c(1,length(cvs))])*(cvs[-c(1,length(cvs))]<cvs[-c(1:2)])==TRUE)+1
   locMaxid <- which((cvs[-c(length(cvs),length(cvs)-1)]< cvs[-c(1,length(cvs))])*(cvs[-c(1,length(cvs))]>cvs[-c(1:2)])==TRUE)+1
   locMinid <- locMinid[locMinid<ifelse(length(locMaxid)>0,max(locMaxid),length(cvs))]
   if(length(locMinid)>0)   opt.lambda=ifelse(  length(locMinid)>2,cand.lambda[ locMinid[which.min(cvs[locMinid[1:(length(locMinid)-1)]])] ],cand.lambda[locMinid[ which.min(cvs[locMinid[1:max(1,length(locMinid))]])] ]  )
   return(opt.lambda/ifelse(length(unique(mscale))==1,4,1.5))
   }


twostep.Binomial <- function(Gramat1,Gramat2,y,wt,lambda0,M)
  {
  n=length(y)
  p=length(wt)

       #--- Step 1.1: All theta=1 ---#
       thetahat=rep(1,p)
       #--- Step 1.2: Solve c,b  ----#
       ssplineGLMobj=sspline.Binomial(Gramat1,Gramat2,thetahat/wt^2,y,lambda0)
       f.old=ssplineGLMobj$fit

       #--- Step 2.1: Solve theta ---#
       pi.old = 1/(1+exp(-f.old))
       W.old = diag(as.numeric(pi.old*(1-pi.old)))
       z = f.old+(y-pi.old)/diag(W.old)
       thetahat=garrote.Logistic(Gramat1,Gramat2,y,wt,ssplineGLMobj$coefs,ssplineGLMobj$intercept,lambda0,M)
       #--- Step 2.2: Solve c,b  again ----#
       ssplineGLMobj=sspline.Binomial(Gramat1,Gramat2,thetahat/wt^2,y,lambda0)
       f.new=ssplineGLMobj$fit

  obj=list(intercept=ssplineGLMobj$intercept,coefs=ssplineGLMobj$coefs,fit=ssplineGLMobj$fit,theta=thetahat)
  return(obj)
  }

garrote.Logistic<- function(Gramat1,Gramat2,y,wt,c0,b0,lambda0,M,init.Theta)
   {
    n=length(y)
    nbasis=dim(Gramat2)[1]
    p=length(wt)
    if(missing(init.Theta))
       {
       L2norm=rep(NA,p)
       for(j in 1:p)  L2norm[j]=sqrt( mean( (wt[j]^(-2)*Gramat1[,,j]%*%c0)^2 ) )
       init.Theta=L2norm*M/sum(L2norm)
       }

     G1=matrix(NA,ncol=p,nrow=n)
     G2=matrix(NA,ncol=p,nrow=nbasis)
     for(j in 1:p)
        {
        G1[,j] =wt[j]^(-2)*Gramat1[,,j]%*%c0
        G2[,j] =wt[j]^(-2)*Gramat2[,,j]%*%c0
        }
     outerG1=array(NA,dim=c(p,p,n))
     for(i in 1:n)  outerG1[,,i]=G1[i,]%*%t(G1[i,])

     Amat = rbind(diag(1,p),rep(-1,p))
     bvec = c(rep(0,p),-M)

     loop=0; Iter.Diff=Inf
     old.theta=init.Theta
     while(loop<10 & Iter.Diff>1e-4)
        {
         GH=garrote.Logistic.GH(G1,G2,y,c0,b0,old.theta,outerG1,lambda0)
         #----- Quadratic Programming -------#
         new.theta <- solve.QP(GH$H,GH$H%*%old.theta-GH$G,t(Amat),bvec)$solution
         new.theta[new.theta<1e-8]<- 0
         #-----------------------------------#
         Iter.Diff=max(abs(new.theta-old.theta))
         old.theta=new.theta
         loop=loop+1
        }
   return(new.theta)
   }

garrote.Logistic.GH=function(G1,G2,y,c0,b0,init.Theta,outerG1,lambda0)
  {
  link=b0+G1%*%init.Theta
  pihat=as.numeric(exp(link)/(1+exp(link)))
  G=apply( diag( pihat-y )%*%G1,2,sum ) +length(y)*lambda0*t(G2)%*%c0
  H=apply( outerG1*array(rep(pihat*(1-pihat),each=length(init.Theta)^2),dim=c(length(init.Theta),length(init.Theta),length(y))),c(1,2),sum )
  return(list(G=G,H=H))
  }

sspline.Binomial <- function(Gramat1,Gramat2,mscale,y,lambda0)
  {
  n=length(y)
  nbasis=dim(Gramat2)[1]
  p=length(mscale)

  Rtheta1=wsGram(Gramat1,mscale)
  Rtheta2=wsGram(Gramat2,mscale)
  EigRtheta2=eigen(Rtheta2);
  loop=0
  while(min(EigRtheta2$values)<0 & loop<10)
     {
     loop=loop+1
     Rtheta2=Rtheta2+1e-8*diag(nbasis)
     EigRtheta2=eigen(Rtheta2);
     }
  if(loop==10)   EigRtheta2$values[EigRtheta2$values<0]=1e-8
  #----- Use glmnet to solve coefficients C -----$
  pseudoX=Rtheta1%*%EigRtheta2$vectors%*%diag(sqrt(1/EigRtheta2$values))
  ssGLM.en=glmnet(pseudoX,y,family="bin",lambda=c(lambda0/2,lambda0,lambda0*2),alpha=0,standardize = FALSE)
  ssGLM.en.coef=as.matrix(predict(ssGLM.en,s=lambda0,type="coef"))
  chat=as.numeric( EigRtheta2$vectors%*%diag(sqrt(1/EigRtheta2$values))%*%ssGLM.en.coef[-1] )
  fhat=Rtheta1%*%chat+ssGLM.en.coef[1]
  return(list(intercept=ssGLM.en.coef[1],coefs=chat,Df=NULL,fit=fhat))
  }

SSANOVAwt.Binomial <- function(x,y,mscale,nbasis,basis.id)
  {
  n <- nrow(x)
  p <- ncol(x)
  if( missing(nbasis) &  missing(basis.id))
    {
    nbasis=max(45, ceiling(12 * length(time)^(2/9)))
    basis.id=sort(sample(1:n,nbasis))
    }
  if( missing(nbasis) & !missing(basis.id))     nbasis <- length(basis.id)
  if(!missing(nbasis) &  missing(basis.id))     basis.id <- sort(sample(1:n,nbasis))

  GramatF <- bigGram(x,x)
  Gramat1 <- GramatF[        ,basis.id,]
  Gramat2 <- GramatF[basis.id,basis.id,]

  optLambda0=cvlam.logistic(Gramat1,Gramat2,y,mscale,10,2^(seq(-8,-21,-1)),T)
  ssObj=sspline.Binomial(Gramat1,Gramat2,mscale,y,optLambda0)
  L2norm=rep(NA,p)
  for(j in 1:p)  L2norm[j]=ifelse(length(unique(x[,j]))>6,sqrt(mean((mscale[j]*Gramat1[,,j]%*%ssObj$coefs)^2)),diff(range(mscale[j]*Gramat1[,,j]%*%ssObj$coefs)))

  return(1/L2norm)
  }

