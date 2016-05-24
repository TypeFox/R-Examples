plot.svmpath=function(x,step=max(x$Step),Size=60,elbow.show=TRUE,support.show=TRUE,...){
  object=x
   if(!missing(step)){
    maxstep=max(object$Step)
    if(step>maxstep){
      warning(paste("maximim step of",format(trunc(maxstep)),"used"))
      step=maxstep
    }
    if(step<=0){
       warning("minimim step of 1 used")
       step=1
     }
  }
  f=predict(object,lambda=object$lambda[step],type="function")
  x=object$x
  y=object$y
  Elbow=object$Elbow[[step]]
  stats<-StatPath(y,f,Elbow)
  alpha=object$alpha[,step]
  alpha0=object$alpha0[step]
  lambda=object$lambda[step]
  linear.plot=object$linear
 ###only for 2 dim inputs
  if(dim(x)[2]>2){
    warning("plot.svm intended for two-dimensional x; only first two dimensions used")
    x <- x[, 1:2]
  }
  n<-length(y)
    ss<-abs(diff(range(x[,2])))
    soffset<-x*0
    soffset[,2]<-ss/50
    x <- x[, 1:2]
  
if(!linear.plot){
    kernel.function=object$kernel
    param.kernel=object$param.kernel
    xr <- apply(x, 2, range)
    xg <- apply(xr, 2, function(x, Size)
                seq(from = x[1], to = x[2], length = Size), Size)
    xG <- data.matrix(expand.grid(as.list(data.frame(xg))))
    Knew <- kernel.function(x, xG, param.kernel=param.kernel,...)
    fhat <- ((alpha * y) %*% Knew + alpha0)/lambda
      fhat <- matrix(fhat, Size, Size)
      }
    dotts=list(...)
    plotargs=list(x=x[,  ], type = "n",...,xlab="X1",ylab="X2",
    main=paste("Step: ",format(step,digits=3)," Errors:",format(round(stats$error),digits=3), " Elbow Size:",format(round(stats$selbow),digits=2)," Sum Eps:",format(round(stats$margin,2),digits=7)))
    plotargs[names(dotts)]=NULL
    plotargs[names(dotts)]=dotts
    do.call("plot",plotargs)
    pointid <- seq(y)
    support=(y*f <= 1)|match(pointid,Elbow,FALSE)
   if(!support.show)    points(x, col = y+3,pch="*",cex=2)
  else{
    points(x[support,], col = y[support]+3,pch=42,cex=2)
    points(x[!support,], col = y[!support]+3,pch=1,cex=1)
  }
    if(n<15)text((x-soffset), labels = paste(pointid), col = y+3)
    if(n<15&&length(Elbow))text((x-soffset)[Elbow,  ], labels = paste(pointid[Elbow]), col = 3)
  if(linear.plot){
    beta <- (alpha * y) %*% x
    abline( - alpha0/beta[2],  - beta[1]/beta[2], col = 3, lwd = 3)
    abline(lambda/beta[2] - alpha0/beta[2],  - beta[1]/beta[2], col = 3, 
           lwd = 2, lty = 3)
    abline( - lambda/beta[2] - alpha0/beta[2],  - beta[1]/beta[2], col = 3, 
           lwd = 2, lty = 3)
  }
  else{
      contour(xg[, 1], xg[, 2], fhat, levels = 0, add = TRUE, labels = "",col=3,lwd=3)
      contour(xg[, 1], xg[, 2], fhat, levels = c(-1,1), add = TRUE, labels = c("",""),col=3,lwd=2,lty=3)
    }
  if(elbow.show)    points(x[Elbow,], col = y[Elbow]+3,pch=20,cex=2.2)

    invisible(list(x=x,y=y,Elbow=Elbow,alpha=alpha,alpha0=alpha0,support=pointid[support],step=step,error=stats$error,sum.eps=stats$margin))
}
