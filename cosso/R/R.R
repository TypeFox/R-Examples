.packageName <- "cosso"

#==================================#
#------- Generic Functions --------#
SSANOVAwt <- function(x,y,tau,family=c("Gaussian","Binomial","Cox","Quantile"),mscale=rep(1,ncol(x)),gamma=1,scale=FALSE,nbasis,basis.id,cpus)
   {
    family = match.arg(family)
    n <- nrow(x)
    p <- ncol(x)
    if(scale) x <- apply(x,2,rescale)
    if(missing(cpus)) cpus=1
    if(family=="Cox" & !all(match(c("time", "status"), dimnames(y)[[2]], 0)) ) 
       {
       stop("Cox model requires a matrix with columns 'time' and 'status' as a response")
       }
    if(family=="Quantile")  
       {
       if((tau<=0)+(tau>=1))  stop("Quantile model requires the value tau strictly between 0 and 1")
       }

   inv.L2norm=switch(family,Gaussian=SSANOVAwt.Gaussian(x,y,mscale,basis.id=basis.id),
                        Cox     =SSANOVAwt.Cox(x,y[,"time"],y[,"status"],mscale),
                        Binomial=SSANOVAwt.Binomial(x,y,mscale,nbasis,basis.id),
                        Quantile=SSANOVAwt.qr(x,y,tau,ifelse(n>=200,10,5),cpus) )
   
   return(inv.L2norm^gamma)
   }

cosso <- function(x,y,tau,family=c("Gaussian","Binomial","Cox","Quantile"),wt=rep(1,ncol(x)),scale=FALSE,nbasis,basis.id,cpus)
   {
    family = match.arg(family)
    n <- nrow(x)
    p <- ncol(x)
    if(scale) x <- apply(x,2,rescale)
    if(missing(cpus)) cpus=1
    if(family=="Cox" & !all(match(c("time", "status"), dimnames(y)[[2]], 0)) ) 
       {
       stop("Cox model requires a matrix with columns 'time' and 'status' as a response")
       }
    if(family=="Quantile")  
       {
       if((tau<=0)+(tau>=1))  stop("Quantile model requires the value tau strictly between 0 and 1")
       }
    GramatF <- bigGram(x,x)

    obj=switch(family,Gaussian=cosso.Gaussian(GramatF,y,wt,nbasis,basis.id),
                      Cox=cosso.Cox(          GramatF,y[,"time"],y[,"status"],wt,nbasis,basis.id,ifelse(cpus>1,TRUE,FALSE),cpus),
                      Binomial=cosso.Binomial(GramatF,y,wt,nbasis,basis.id),
                      Quantile=cosso.qr(      GramatF,y,tau,wt,ifelse(cpus>1,TRUE,FALSE),cpus)  )
    obj=c(list(y=y,x=x,Kmat=GramatF,wt=wt),obj)                  
    return(obj)
   }
   
tune.cosso <- function(object,folds=5,plot.it=TRUE)   
   {
   family=object$family
   if (folds <= 3)       stop("folds must be bigger than or equal to 3; folds=5 or 10 recommended")   

   tuningobj=switch(family,Gaussian=tune.cosso.Gaussian(object,folds,plot.it),
                           Cox     =tune.cosso.Cox(     object,folds,plot.it),
                           Binomial=tune.cosso.Binomial(object,folds,plot.it),
                           Quantile=tune.cosso.qr(      object,folds,plot.it))
   return(tuningobj)                  
   }

plot.cosso<- function(x,M,plottype =c("Path","Functionals"),eps=1e-7,...)
  {
  object <- x
  p<- ncol(object$x)
  n<- nrow(object$x)
  plottype <- match.arg(plottype)
  if (missing(M) & plottype == "Functionals")
     {
      warning("Type=fit with no smoothing parameter argument; type switched to Path")
      plottype <- "Path"
     }
     
  if(plottype=="Path")
     {
      matplot(object$tune$Mgrid,object$tune$L2norm,type="l",lty=1,col=c(1,rainbow(p-1)),xlab="M",ylab=expression(L[2]-norm))
      axis(4,at=object$tune$L2norm[length(object$tune$Mgrid),],labels=1:p,cex=.6,las=2)
     }
  else
     {
      fitObj=switch(object$family,Gaussian=twostep.Gaussian(object$Kmat[,object$basis.id,],object$Kmat[object$basis.id,object$basis.id,],object$y,object$wt,object$tune$OptLam,M),
                                  Cox     =twostep.Cox(     object$Kmat[,object$basis.id,],object$Kmat[object$basis.id,object$basis.id,],object$y[,"time"],object$y[,"status"],object$wt,object$RiskSet,object$tune$OptLam,M),
                                  Binomial=twostep.Binomial(object$Kmat[,object$basis.id,],object$Kmat[object$basis.id,object$basis.id,],object$y,object$wt,object$tune$OptLam,M),
                                  Quantile=twostep.qr(object$tau,object$y,object$Kmat,object$tune$OptLam,M,object$wt)  )

      selectID <- (1:p)[fitObj$theta>eps]
      funcProfile=matrix(NA,ncol=p,nrow=n)
      for (j in selectID)   funcProfile[,j]=fitObj$theta[j]/(object$wt[j]^2)*object$Kmat[,object$basis.id,j]%*%fitObj$coefs

      plotrow <- ceiling(length(selectID)/2)
      plotcol <- 2

      par(mfrow=c(plotrow,plotcol))
      for (j in selectID) 
        {
        if(is.null(colnames(object$x)[j]))  xlab=paste("x",j,sep="")
        else                                xlab=colnames(object$x)[j]
        ylab=paste("f(",xlab,")",sep="")
        categorical=length(unique(object$x[,j]))<7

        if(categorical)
            {
            plotData= unique(cbind(object$x[,j],funcProfile[,j]))
            plot(0,0,col=0,xlab=xlab,ylab=ylab,xlim=range(plotData[,1])+c(-0.4,0.4),ylim=1.1*range(plotData[,2]))
            for(k in 1:nrow(plotData))  segments(plotData[k,1]-0.3,plotData[k,2],plotData[k,1]+0.3,plotData[k,2],lwd=2)
            }
        else
            {
            xorder=order(object$x[,j])
            plot(object$x[xorder,j],funcProfile[xorder,j],type="l",xlim=c(0,1),xlab=xlab,ylab=ylab)
            }
         }
     }
  }


predict.cosso<- function(object,xnew,M,type=c("fit","coefficients","nonzero"),eps=1e-7,...)
  {
  type<- match.arg(type)
  if (missing(xnew) & type == "fit")
     {
      warning("Type=fit with no xnew argument; type switched to coefficients")
      type <- "coefficients"
     }
  if(missing(M))
    {
     warning("Missing Smoothing Parameter, M. Replaced by value tuned by 5-fold CV")
     tuneObj <- tune.cosso(object,folds=5,plot.it=FALSE) 
     M <- tuneObj$OptM
    }

  fitObj=switch(object$family,Gaussian=twostep.Gaussian(object$Kmat[,object$basis.id,],object$Kmat[object$basis.id,object$basis.id,],object$y,object$wt,object$tune$OptLam,M),
                              Cox     =twostep.Cox(     object$Kmat[,object$basis.id,],object$Kmat[object$basis.id,object$basis.id,],object$y[,"time"],object$y[,"status"],object$wt,object$RiskSet,object$tune$OptLam,M),
                              Binomial=twostep.Binomial(object$Kmat[,object$basis.id,],object$Kmat[object$basis.id,object$basis.id,],object$y,object$wt,object$tune$OptLam,M),
                              Quantile=twostep.qr(object$tau,object$y,object$Kmat,object$tune$OptLam,M,object$wt)  )

  if (type=="fit")
     {
      predictor = as.numeric(fitObj$intercept+wsGram(bigGram(xnew, object$x[object$basis.id,]), fitObj$theta/(object$wt^2))%*%fitObj$coefs)
     }
  else if(type=="coefficients")
     {
      predictor = list(coef=fitObj$coefs,intercept=fitObj$intercept,theta=fitObj$theta)
     }
  else
     {
     predictor=which(fitObj$theta>eps)
     }
  class(predictor)="predict.cosso"
  return(predictor)
  }

