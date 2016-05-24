plot.ITP1 <-
function(x,xrange=c(0,1),alpha1=0.05,alpha2=0.01,
                      ylab='Functional Data',main=NULL,lwd=1,col=1,pch=16,ylim=range(object$data.eval),
                      ...){
  if(alpha1 < alpha2){
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
  object <- x
  
  par(ask=TRUE) 
  if(object$basis=='Fourier'){
    p <- length(object$pval)
    J <- dim(object$data.eval)[2]
    n <- dim(object$data.eval)[1]
    xmin <- xrange[1]
    xmax <- xrange[2]
    abscissa.pval = 1:p
    Abscissa = seq(xmin,xmax,len=J)
    main.data <- paste(main,': Functional Data')
    main.data <- sub("^ : +", "", main.data)
    matplot(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=col,lwd=lwd,ylim=ylim,...)
    if(length(object$mu)==1){
      abscissa.mu <- Abscissa
      mu <- rep(object$mu,1000)
    }else{
      Abscissa <- seq(xmin,xmax,length.out=length(object$mu))
      mu <- object$mu
    }
    lines(abscissa.mu,mu,col='gray',lwd=2)
    
    ################################################################
    # pval
    main.p <- paste(main,': Adjusted p-values')
    main.p <- sub("^ : +", "", main.p)
    plot(abscissa.pval,object$corrected.pval,pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',xlab='Frequency')
    difference1 <- which(object$corrected.pval<alpha1)
    if(length(difference1)>0){
      for(j in 1:length(difference1)){
        min.rect <- abscissa.pval[difference1[j]] - 0.5
        max.rect <- min.rect + 1
        rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = 'gray90',density=-2,border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = NULL,border='black')
    }
    difference2 <- which(object$corrected.pval<alpha2)
    if(length(difference2)>0){
      for(j in 1:length(difference2)){
        min.rect <- abscissa.pval[difference2[j]] - 0.5
        max.rect <- min.rect + 1
        rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = 'gray80',density=-2,border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = NULL,border='black')
    }
    for(j in 0:10){
      abline(h=j/10,col='lightgray',lty="dotted")
    }
    points(1:p,object$corrected.pval,pch=pch)
    
    
    
  }else if(object$basis=='B-spline'){  
    p <- length(object$pval)
    J <- dim(object$data.eval)[2]
    n <- dim(object$data.eval)[1]
    xmin <- xrange[1]
    xmax <- xrange[2]
    abscissa.pval = seq(xmin,xmax,len=p)
    Abscissa = seq(xmin,xmax,len=J)
    main.data <- paste(main,': Functional Data')
    main.data <- sub("^ : +", "", main.data)
    matplot(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=col,lwd=lwd,ylim=ylim,...)
    difference1 <- which(object$corrected.pval<alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray90", density = -2, border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
    }
    difference2 <- which(object$corrected.pval<alpha2)
    if (length(difference2) > 0) {
      for (j in 1:length(difference2)) {
        min.rect <- abscissa.pval[difference2[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray80", density = -2, border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
    }
    matplot(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=col,lwd=lwd,add=TRUE,...)
    if(length(object$mu)==1){
      abscissa.mu <- Abscissa
      mu <- rep(object$mu,1000)
    }else{
      Abscissa <- seq(xmin,xmax,length.out=length(object$mu))
      mu <- object$mu
    }
    lines(abscissa.mu,mu,col='gray',lwd=2)
    
    ################################################################
    # pval
    main.p <- paste(main,': Adjusted p-values')
    main.p <- sub("^ : +", "", main.p)
    plot(abscissa.pval,object$corrected.pval,pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',...)
    difference1 <- which(object$corrected.pval<alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray90", density = -2, border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
    }
    difference2 <- which(object$corrected.pval<alpha2)
    if (length(difference2) > 0) {
      for (j in 1:length(difference2)) {
        min.rect <- abscissa.pval[difference2[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray80", density = -2, border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
    }
    for(j in 0:10){
      abline(h=j/10,col='lightgray',lty="dotted")
    }
    points(abscissa.pval,object$corrected.pval,pch=pch)
    
  }
  par(ask=FALSE) 
}
