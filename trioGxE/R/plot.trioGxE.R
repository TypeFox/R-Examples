plot.trioGxE <- function(x,se=TRUE,seWithGxE.only=TRUE,ylim=NULL,
                         yscale=TRUE,xlab=NULL,ylab=NULL,rugplot=TRUE,...){
  ## ----------------------------------------------------------------------
  ## x: returned object from pirls.trio()
  ## se: whether to add the estimated variance for the curve estimates.
  ##
  ## seWithGxE.only = If TRUE, the estimated intervals for GxE curves account
  ## for uncertainty in the centered smoothed GattrE curves only. If FALSE, the 
  ## Bayesian credible intervals include the uncertainty about the main genetic 
  ## effect as well.
  ##
  ## yscale=plot the two GRR curves on the same y-scale??
  ## 
  ## ...: arguments passed to par() function 
  ## ----------------------------------------------------------------------
  mt <- x$triodata$mt; ind.sm1<- mt==1|mt==3.1; ind.sm2<- mt==2|mt==3.2
  attr <- x$triodata$x #augmented (mt1+mt2+mt3.1+mt3.2) covariate vector
  x1 <- attr[ind.sm1]; x2 <- attr[ind.sm2]
  
  ## calculating confidence intervals
  smooth.ci(object=x,se=se,seWithGxE.only=seWithGxE.only,
            ylim=ylim,yscale=yscale,xlab=xlab,ylab=ylab,mfrow=c(1,2),
            rugplot=rugplot,x1=x1,x2=x2,ind.sm1=ind.sm1,ind.sm2=ind.sm2,...)
  
}

smooth.ci <- function(object,se,seWithGxE.only,ylim,yscale,
                      xlab,ylab,mfrow,rugplot,
                      x1,x2,ind.sm1,ind.sm2,...){
  ## ... includes parameter that should be passed to par() ??
  ## lwd=1,cex.main=1.5,cex.lab=1.25,mar=c(4,5,3,1)
  ## extracting variance-covariance matrix from the pirls.trio object
  Var.coef = object$Vp
  
  ## extract basis info
  penmod = object$penmod
  k1 <- object$smooth$bs.dim[1]
  k2 <- object$smooth$bs.dim[2]
  xk1 <- object$smooth$knots[[1]]
  xk2 <- object$smooth$knots[[2]]
  
  ## set the confidence margin: default is 2 se above and below the estimated curves
  if(se > 0 ){ #se is either TRUE or a positive number
    if( se == TRUE ) #default, 2 se curves are plotted
      err.mult = 2   
    else # a positive number
      err.mult = se 
  }
  
  else #se is FALSE or a negative number
    err.mult = 0

  range.x1 = range(x1); range.x2 = range(x2)
  
  ## Caculating GRRs for the new set of grid-points 
  if(penmod =="codominant"){
    
    tem.x1 = seq(from=range.x1[1],to=range.x1[2],length.out=100)
    tem.x2 = seq(from=range.x2[1],to=range.x2[2],length.out=100)
    
    ## GRR_1(e): genotype risk for G=1 vs. G=0
    X1 = Xmat(xK=xk1, x=tem.x1) 
    qrc <- object$qrc$qrc1
    if(!is.null(qrc)){
      tem.X <- X1
      tem.X[,-1] <- t(qr.qy(qrc,t(X1))[2:k1,])
      tem.X[,1] <- 1
      X1 <- tem.X
      rm(tem.X)
    }
    else{
      tem.X <- cbind(1,X1)
      X1 <- tem.X
      rm(tem.X)
    }
    rm(qrc)
    coef1 <- object$coef[1:k1] 
    Var.coef1=Var.coef[1:k1,1:k1]     
    
    ## GRR_2(e): genotype risk for G=2 vs. G=1
    X2 = Xmat(xK=xk2, x=tem.x2) 
    qrc <- object$qrc$qrc2
    if(!is.null(qrc)){
      tem.X <- X2
      tem.X[,-1] <- t(qr.qy(qrc,t(tem.X))[2:k2,])
      tem.X[,1] <- 1
      X2 <- tem.X
      rm(tem.X)
    }
    else{
      tem.X <- cbind(1,X2)
      X2 <- tem.X
      rm(tem.X)
    }
    rm(qrc)

    coef2 <- object$coef[-c(1:k1)] #basis coefficient estimates
    Var.coef2=Var.coef[-c(1:k1),-c(1:k1)]
    
    ## plot without genetic main effect terms
    if(seWithGxE.only){
      Var.curve1 = X1[,-1]%*%Var.coef1[-1,-1]%*%t(X1[,-1])
      Var.curve2 = X2[,-1]%*%Var.coef2[-1,-1]%*%t(X2[,-1])    
    }
    
    else{
      Var.curve1 = X1%*%Var.coef1%*%t(X1)
      Var.curve2 = X2%*%Var.coef2%*%t(X2)
    }
    
    fit.curve1 <- X1%*%coef1 - coef1[1]#mean.f1 ##centered (i.e., f1.hat)
    fit.curve2 <- X2%*%coef2 - coef2[1]#mean.f2 ##centered (i.e., f2.hat)
    
    ## calculate upper and lower confidence limits for a given margin
    ci.up1 <- fit.curve1 + rep(err.mult,100)*sqrt(diag(Var.curve1))
    ci.up2 <- fit.curve2 + rep(err.mult,100)*sqrt(diag(Var.curve2))
    
    ci.low1 <- fit.curve1 - err.mult*sqrt(diag(Var.curve1))
    ci.low2 <- fit.curve2 - err.mult*sqrt(diag(Var.curve2))
    
    plot.pirls.gam(object=object,se=se,seWithGxE.only=seWithGxE.only,ylim=ylim,yscale=yscale,
                   xlab=xlab,ylab=ylab,mfrow=mfrow,rugplot=rugplot,
                   tem.x1=tem.x1,tem.x2=tem.x2,x1=x1,x2=x2,k1=k1,k2=k2,
                   fit.curve1=fit.curve1,ci.up1=ci.up1,ci.low1=ci.low1,
                   fit.curve2=fit.curve2,ci.up2=ci.up2,ci.low2=ci.low2,...)
  }#if: codominant penmod

  else{ #non-codominant penmod
    if(penmod=="dominant"){

      tem.x1 = seq(from=range.x1[1],to=range.x1[2],length.out=100)
      tem.x2 = seq(from=range.x2[1],to=range.x2[2],length.out=100)
      
      X1 = Xmat(xK=xk1, x=tem.x1) ## new data matrix
      tem.X = X1
      qrc <- object$qrc$qrc1
      ## when qrc is not null
      X1[,-1] <- t(qr.qy(qrc,t(tem.X))[2:k1,])
      X1[,1] <- 1
      rm(tem.X,qrc)
      
      coef1 <- object$coef[1:k1] #basis coefficient estimates
      Var.coef1=Var.coef[1:k1,1:k1] #** debugging
      
      ## variance-covariance matrix
      if(seWithGxE.only){
        Var.curve1 = X1[,-1]%*%Var.coef1[-1,-1]%*%t(X1[,-1])
        Var.curve2 = NULL
      }
      else{
        Var.curve1 = X1%*%Var.coef1%*%t(X1)
        Var.curve2 = NULL
      }
      ## fitted values at the new (plotting grid) points
      fit.curve1 <- X1%*%coef1- coef1[1]#mean.f1 ##centered
      fit.curve2 <- rep(0,length(tem.x2)) #object$fitted.values[[2]] #zero
      
      ## calculate the upper and lower confidence limits
      ci.up1 <- fit.curve1 + rep(err.mult,100)*sqrt(diag(Var.curve1))
      ci.up2 <- NULL
    
      ci.low1 <- fit.curve1 - err.mult*sqrt(diag(Var.curve1))
      ci.low2 <- NULL
    }#pen==dom
    
    else if(penmod=="additive"){
      tem.x1 = seq(from=range.x1[1],to=range.x1[2],length.out=100)
      tem.x2 = seq(from=range.x2[1],to=range.x2[2],length.out=100)
      
      X1 = Xmat(xK=xk1, x=tem.x1) ## new data matrix
      qrc <- object$qrc
      if(!is.null(qrc)){
        tem.X <- X1
        tem.X[,-1] <- t(qr.qy(qrc,t(X1))[2:k1,])
        tem.X[,1] <- 1
        X1 <- tem.X
        rm(tem.X)
      }
      else{
        tem.X <- cbind(1,X1)
        X1 <- tem.X
        rm(tem.X)
      }
      
      X2 = Xmat(xK=xk2, x=tem.x2) ## new data matrix
      if(!is.null(qrc)){
        tem.X <- X2
        tem.X[,-1] <- t(qr.qy(qrc,t(tem.X))[2:k2,])
        tem.X[,1] <- 1
        X2 <- tem.X
        rm(tem.X)
      }
      else{
        tem.X <- cbind(1,X2)
        X2 <- tem.X
        rm(tem.X)
      }
      rm(qrc)
      coef1 <- coef2 <- object$coef[1:k1] #basis coefficient estimates
      Var.coef1 <- Var.coef2 <- Var.coef[1:k1,1:k1] #** debugging
      
      if(seWithGxE.only){
        Var.curve1 = X1[,-1]%*%Var.coef1[-1,-1]%*%t(X1[,-1])
        Var.curve2 = X2[,-1]%*%Var.coef2[-1,-1]%*%t(X2[,-1])
      }
      
      else{
        Var.curve1 = X1%*%Var.coef1%*%t(X1)
        Var.curve2 = X2%*%Var.coef2%*%t(X2)
      }
      
      fit.curve1 <- X1%*%coef1 - coef1[1]#mean.f1 ##centered (i.e., f1.hat)
      fit.curve2 <- X2%*%coef2 - coef2[1]#mean.f2 ##centered (i.e., f2.hat)
      
      ci.up1 <- fit.curve1 + rep(err.mult,100)*sqrt(diag(Var.curve1))
      ci.up2 <- fit.curve2 + rep(err.mult,100)*sqrt(diag(Var.curve2))
      
      ci.low1 <- fit.curve1 - err.mult*sqrt(diag(Var.curve1))
      ci.low2 <- fit.curve2 - err.mult*sqrt(diag(Var.curve2))
      
    }#additive
    
    else{#recessive      
      tem.x1 = seq(from=range.x1[1],to=range.x1[2],length.out=100)
      tem.x2 = seq(from=range.x2[1],to=range.x2[2],length.out=100)
      
      X2 = Xmat(xK=xk2, x=tem.x2) ## new data matrix
      tem.X = X2
      qrc <- object$qrc$qrc2
      ## when qrc is not null
      X2[,-1] <- t(qr.qy(qrc,t(tem.X))[2:k2,])
      X2[,1] <- 1
      rm(tem.X,qrc)
      
      coef2 <- object$coef[1:k2] #basis coefficient estimates
      Var.coef2=Var.coef[1:k2,1:k2] #** debugging
      
      ## variance-covariance matrix
      if(seWithGxE.only){
        Var.curve1 = NULL
        Var.curve2 = X2[,-1]%*%Var.coef2[-1,-1]%*%t(X2[,-1])
      }
      else{
        Var.curve1 = NULL
        Var.curve2 = X2%*%Var.coef2%*%t(X2)
      }
      ## fitted values at the new (plotting grid) points
      fit.curve1 <- rep(0,length(tem.x1)) #object$fitted.values[[1]] 
      fit.curve2 <- X2%*%coef2- coef2[1]#mean.f1 ##centered
      
      ci.up1 <- NULL
      ci.up2 <- fit.curve2 + rep(err.mult,100)*sqrt(diag(Var.curve2))
      
      ci.low1 <- NULL
      ci.low2 <- fit.curve2 - err.mult*sqrt(diag(Var.curve2))      
    }#i.e.,(penmod=="recessive")
    ##plotting
    plot.pirls.gam(object=object,se=se,seWithGxE.only=seWithGxE.only,ylim=ylim,yscale=yscale,
                   xlab=xlab,ylab=ylab,mfrow=mfrow,rugplot=rugplot,
                   tem.x1=tem.x1,tem.x2=tem.x2,x1=x1,x2=x2,k1=k1,k2=k2,
                   fit.curve1=fit.curve1,ci.up1=ci.up1,ci.low1=ci.low1,
                   fit.curve2=fit.curve2,ci.up2=ci.up2,ci.low2=ci.low2,...)
  }#!is.null(penmod)
} #smooth.ci() ends here

plot.pirls.gam <- function(object,se,seWithGxE.only,ylim,yscale,xlab,ylab,mfrow,rugplot,
                           fit.curve1,ci.up1,ci.low1,fit.curve2,ci.up2,ci.low2,
                           tem.x1=tem.x1,tem.x2=tem.x2,x1,x2,k1,k2,...)
{
  ## used in smooth.ci() function
  ##object = pirls object?
  ## lwd=1,cex.main=1.5,cex.lab=1.25,mar=c(4,5,3,1)
  
  xlim <- range(c(x1,x2))
  
  ##setting y-limits
  if(is.null(ylim)){   
    if(!se){ # confidence bands not plotted:
      if(yscale){
        ylim1=ylim2=range(c(fit.curve1,fit.curve2))
      }   
      else
        {
          ylim1=range(fit.curve1)
          ylim2=range(fit.curve2)
        }
    }# if(!se)
    else{ # se=TRUE: confidence bands are plotted
      if(yscale){
        ylim1=ylim2=range(c(ci.up1,ci.up2,ci.low1,ci.low2))
      }
      else{ #y-scale is not the same for the two GxE curves
        ylim1=range(c(ci.low1,ci.up1))
        ylim2=range(c(ci.low2,ci.up2))
      }
    } #else(!se)
  }# if(is.null(ylim))

  else{# ylim is provided
    if(is.list(ylim)){
      ylim1 = ylim[[1]]
      ylim2 = ylim[[2]]
    }
    else {# vector
      ylim1 = ylim2 = ylim
    }
  }#else to if(is.null(ylim))
  
  if(object$penmod == "codominant")
    edfVal = c(round(sum(object$edf[c(2:k1)]),2),round(sum(object$edf[-c(1:(k1+1))]),2))
  else if(object$penmod == "dominant")
    edfVal = c(round(sum(object$edf[-1]),2),NA) 
  else if(object$penmod == "additive")
    edfVal = rep(round(sum(object$edf[-1]),2),2) 
  else #recessive
    edfVal = c(NA,round(sum(object$edf[-1]),2))

  if(is.null(ylab)){
    ylab1 <- bquote(hat(f[1])~(edf == .(edfVal[1])))
    ylab2 <- bquote(hat(f[2])~(edf == .(edfVal[2])))    
  }
  else{
    ylab1 <- ylab[[1]]
    ylab2 <- ylab[[2]]
  }
  
  if(is.null(mfrow))
    par(mfrow=c(1,2), mar=c(4,5,3,1),...)

  else{
    if(all(mfrow==c(1,2)))
      par(mfrow=mfrow, mar=c(4,5,3,1),...)
    else
      par(mfrow=mfrow,...)      
  }
  
  ##------------ smooth1 ------------##
  plot(tem.x1,fit.curve1,xlim=xlim,ylim=ylim1,type="l",
       xlab="",ylab=ylab1,main="",...)
  title(main=paste("G=1 vs. G=0 (n=",length(x1),")",sep=""))
  
  if(se)
    {
      lines(tem.x1,ci.up1,lty=2,...)
      lines(tem.x1,ci.low1,lty=2,...)
    }
  
  if(rugplot)
    rug(jitter(x1, amount=0.02))
  
  ##------------ smooth2 ------------##
  plot(tem.x2,fit.curve2,xlim=xlim,ylim=ylim2,type="l",
       xlab="",ylab=ylab2,main="",...)
  title(main=paste("G=2 vs. G=1 (n=",length(x2),")",sep="")) 
  if(rugplot)
    rug(jitter(x2, amount=0.01))
  
  if(se)
    {
      lines(tem.x2,ci.up2,lty=2,...)
      lines(tem.x2,ci.low2,lty=2,...)
    }
  
  ## xlab
  if(is.null(xlab)) 
    xlab=paste(object$terms$cenv)
  else
    xlab=xlab

  mtext(xlab,outer=TRUE,cex=par()$cex.lab,side=1,line=-1.25)
  
}#function ends
