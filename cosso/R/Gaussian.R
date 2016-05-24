#######=============================#####
#--- Functions for Gaussian Response ---#
#######=============================#####

cosso.Gaussian <- function(Gramat,y,wt,nbasis,basis.id)
   {
     n <- length(y)
     p <- length(wt)
     if( missing(nbasis) &  missing(basis.id))
       {
        nbasis=max(40, ceiling(12 * n^(2/9)))
        basis.id=sort(sample(1:n,nbasis))
       }
     if( missing(nbasis) & !missing(basis.id))     nbasis   <- length(basis.id)
     if(!missing(nbasis) &  missing(basis.id))     basis.id <- sort(sample(1:n,nbasis))
     Gramat1 <- Gramat[        ,basis.id, ]
     Gramat2 <- Gramat[basis.id,basis.id, ]

     bestlam <- cvlam.Gaussian(Gramat1,Gramat2,y,1/wt^2,folds=6)
     L2normMat <- Mgrid <- NULL
     tempM<- 0.2
     tempTheta=tempL2norm=rep(0,p)
     loop=0
     while( sum(tempTheta>1e-7)<p  &  loop<=ifelse(p<=15,floor(2*p),p) )
         {
         loop=loop+1
         Mgrid  <- c(Mgrid,tempM)
         tmpSol <- twostep.Gaussian(Gramat1,Gramat2,y,wt,bestlam,tempM)
         #--- Compute Solution Path ---#
         for(j in 1:p)  L2normMat<- c( L2normMat,sqrt(mean((tmpSol$theta[j]/wt[j]^2*Gramat1[,,j]%*%tmpSol$coefs)^2)) )
         for(j in 1:p)  tempL2norm[j] <- ifelse( length(unique(as.numeric(Gramat[,,j])))>6 , sqrt(mean((tmpSol$theta[j]/wt[j]^2*Gramat1[,,j]%*%tmpSol$coefs)^2)) , diff(range(tmpSol$theta[j]/wt[j]^2*Gramat1[,,j]%*%tmpSol$coefs)) )

         if(loop<10)                   tempM <- tempM+0.25
         else if(10<=loop & loop<16)   tempM <- tempM+0.5
         else if(16<=loop & loop<20)   tempM <- tempM+1
         else                          tempM <- tempM+2
         }
     L2normMat <- matrix(L2normMat,ncol=p,byrow=TRUE)
     L2normMat <- rbind(rep(0,p),L2normMat)
     Mgrid <- c(0,Mgrid)
     cossoobj<- list(family="Gaussian",basis.id=basis.id,tune=list(OptLam=bestlam,Mgrid=Mgrid,L2norm=L2normMat) )
     class(cossoobj)="cosso"
     return(cossoobj)
   }

tune.cosso.Gaussian <- function(object,folds=5,plot.it=TRUE)
   {
    n <- length(object$y)
    p <- length(object$wt)

    origMgrid  <- object$tune$Mgrid[-1]
    uniqueSize <- unique(apply(object$tune$L2norm[-1,]>0,1,sum))
    newGrid <- origMgrid[apply(object$tune$L2norm[-1,]>0,1,sum)<=ifelse(p<=15,floor(p/2),floor(p/3))]
    for(k in 1:length(uniqueSize))   
       {
       if(uniqueSize[k]>ifelse(p<=15,floor(p/2),floor(p/3)))  newGrid=c(newGrid,  origMgrid[max(which(apply(object$tune$L2norm[-1,]>0,1,sum)==uniqueSize[k]))]  )
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

    IDmat<- cvsplitID(n,folds)
    cand.M <- sort(c(newGrid[-1],extMgrid))
    cvRaw  <- matrix(NA,ncol=length(cand.M),nrow=n)

    for(f in 1:folds)
       {
        testID=IDmat[!is.na(IDmat[,f]),f]
        trainID=(1:n)[-testID]
        trainGramat1=object$Kmat[trainID,object$basis.id,]
        testGramat1 =object$Kmat[ testID,object$basis.id,]
        for(m in 1:length(cand.M))
            {
             tempSol =twostep.Gaussian(trainGramat1,object$Kmat[object$basis.id,object$basis.id,],object$y[trainID],object$wt,object$tune$OptLam,cand.M[m])
             tempfpred = tempSol$intercept+wsGram(testGramat1,tempSol$theta/object$wt^2 )%*%tempSol$coefs
             cvRaw[testID,m]=(object$y[testID]-tempfpred)^2
            }
       }
   cvm =apply(cvRaw,2,mean)
   cvsd=sqrt(apply(scale(cvRaw, cvm, FALSE)^2, 2,mean)/n)

   locMinid <- which((cvm[-c(length(cvm),length(cvm)-1)]> cvm[-c(1,length(cvm))])*(cvm[-c(1,length(cvm))]<cvm[-c(1:2)])==TRUE)+1
   locMaxid <- which((cvm[-c(length(cvm),length(cvm)-1)]< cvm[-c(1,length(cvm))])*(cvm[-c(1,length(cvm))]>cvm[-c(1:2)])==TRUE)+1
   locMaxid <- locMaxid[locMaxid>(length(cand.M)*2/3)]
   locMinid <- locMinid[locMinid<ifelse(length(locMaxid)>0,max(locMaxid),length(cvm))]
   opt.M=cand.M[which.min(cvm)]
   if(length(locMinid)>0)   opt.M=ifelse( length(locMinid)>2,cand.M[ locMinid[which.min(cvm[locMinid[1:(length(locMinid)-1)]])] ],cand.M[locMinid[ which.min(cvm[locMinid[1:max(1,length(locMinid))]])] ] )

   if(plot.it)
       {
        par(mfcol=c(1,2))
        plot(cand.M,cvm,ylim=c(min(cvm-cvsd),max(cvm+cvsd)),type="b",col=2,pch=16,lwd=1.5,xlab="M",ylab="Cross-Validated Squared Error")
        for(m in 1:length(cand.M))  segments(cand.M[m],cvm[m]-cvsd[m],cand.M[m],cvm[m]+cvsd[m],col=grey(0.6))
        abline(v=opt.M,lty=2,col=2);axis(3,opt.M)
        matplot(object$tune$Mgrid,object$tune$L2norm,type="l",lty=1,col=c(1,rainbow(p-1)),xlab="M",ylab=expression(L[2]-norm))
        abline(v=opt.M,lty=2,col=2);axis(3,opt.M)
        axis(4,at=object$tune$L2norm[length(object$tune$Mgrid),],labels=1:p,cex=.3,las=2)
        }
    return(list(OptM=opt.M,M=cand.M,cvm=cvm,cvsd=cvsd))
    }

twostep.Gaussian <- function(Gramat1,Gramat2,y,wt,lambda,mm)
   {
    n <- length(y)
    nbasis <- dim(Gramat1)[2]
    d <- length(wt)
    #---- Step 1.2 ----#
    cb0<- sspline(Gramat1,Gramat2,y,1/(wt^2),lambda)
    c0 <- cb0[-1]
    b0 <- cb0[ 1]
    #---- Step 2.1 ----#
    G1 <- matrix(0,n,d)
    G2 <- matrix(0,nbasis,d)
    for(j in 1:d)
       {
       G1[,j] = Gramat1[,,j]%*%c0*(wt[j]^(-2))
       G2[,j] = Gramat2[,,j]%*%c0*(wt[j]^(-2))
       }
    dvec <- 2*t(G1)%*%(y-b0)-n*lambda*t(G2)%*%c0
    Dmat <- 2*t(G1)%*%G1
    Amat <- rbind(diag(d),rep(-1,d))
    bvec   <- c(rep(0,d),-mm)
    theta1 <- solve.QP(Dmat,dvec,t(Amat),bvec)[[1]]
    theta1[theta1<1e-8] <- 0
    #---- Step 2.2 ----#
    cb1 <- sspline(Gramat1,Gramat2,y,theta1/(wt^2),lambda)

    result <- list(coefs=cb1[-1],intercept=cb1[1],theta=theta1)
    return(result)
  }

sspline <- function(Gramat1,Gramat2,y,mscale,lambda)
  {
    n <- length(y)
    d <- length(mscale)
    nbasis <- dim(Gramat1)[2]

    Rtheta1 = wsGram(Gramat1,mscale)
    Rtheta2 = wsGram(Gramat2,mscale)

    LHS=t(Rtheta1)%*%scale(Rtheta1,scale=F)+2*n*lambda*Rtheta2
    RHS=t(Rtheta1)%*%(y-mean(y))
    chat=solve.singular(LHS,RHS)
    bhat=mean(y-Rtheta1%*%chat)
    return(c(bhat,chat))
 }


cvlam.Gaussian <- function(Gramat1,Gramat2,y,mscale,folds=5,cand.lambda)
  { 
    if(missing(cand.lambda))  cand.lambda <- 2^seq(-10,-22,-.75)
    n<-length(y)
    Rtheta1 <- wsGram(Gramat1,mscale)
    Rtheta2 <- wsGram(Gramat2,mscale)


    IDmat<- cvsplitID(n,folds)
    cvRaw<- matrix(NA,ncol=length(cand.lambda),nrow=n)
    for(f in 1:folds)
        {
        testID  <- IDmat[!is.na(IDmat[,f]),f]
        trainID <- (1:n)[-testID]
        trainRtheta1 <- Rtheta1[trainID,]
        testRtheta1  <- Rtheta1[testID,] 
        LHS <- t(trainRtheta1)%*%trainRtheta1
        RHS <- t(trainRtheta1)%*%(y[trainID]-mean(y[trainID]))
        for(k in 1:length(cand.lambda))
           {
           chat  <- solve.singular(LHS+2*length(trainID)*cand.lambda[k]*Rtheta2,RHS)
           bhat  <- mean(y[trainID]-trainRtheta1%*%chat)
           fpred <-bhat+testRtheta1%*%chat
           cvRaw[testID,k] <- (y[testID]-fpred)^2
           }
       }
   cvm <- apply(cvRaw,2,mean)
   locMinid <- which((cvm[-c(length(cvm),length(cvm)-1)]> cvm[-c(1,length(cvm))])*(cvm[-c(1,length(cvm))]<cvm[-c(1:2)])==TRUE)+1
   locMaxid <- which((cvm[-c(length(cvm),length(cvm)-1)]< cvm[-c(1,length(cvm))])*(cvm[-c(1,length(cvm))]>cvm[-c(1:2)])==TRUE)+1
   locMaxid <- locMaxid[locMaxid>(length(cvm)*2/3)]
   locMinid <- locMinid[locMinid<ifelse(length(locMaxid)>0,max(locMaxid),length(cvm))]
   optLambda=cand.lambda[which.min(cvm)]
   if(length(locMinid)>0)   optLambda=ifelse( length(locMinid)>2,cand.lambda[ locMinid[which.min(cvm[locMinid[1:(length(locMinid)-1)]])] ],cand.lambda[locMinid[ which.min(cvm[locMinid[1:max(1,length(locMinid))]])] ] )

   return(optLambda)
  }



SSANOVAwt.Gaussian <- function(x,y,mscale,nbasis,basis.id)
  { 
    n <- length(y)
    d <- ncol(x)
    if( missing(nbasis) &  missing(basis.id))
       {
        nbasis=max(40, ceiling(12 * n^(2/9)))
        basis.id=sort(sample(1:n,nbasis))
       }
    else if(!missing(nbasis) &  missing(basis.id))     basis.id <- sort(sample(1:n,nbasis))
    Gramat1 <- bigGram(x,x[basis.id,])
    Gramat2 <- Gramat1[basis.id,,]

    bestlam <- cvlam.Gaussian(Gramat1,Gramat2,y,mscale,folds=8)
    cc <- sspline(Gramat1,Gramat2,y,mscale,bestlam)[-1]
    L2norm <- rep(NA,d)
    for (j in 1:d)
      {
      fjj=mscale[j]*Gramat1[,,j]%*%cc
      L2norm[j] <- ifelse(length(unique(fjj))>6,sqrt(mean(fjj^2)),diff(range(fjj)))
      }
    return(1/L2norm)
  }

