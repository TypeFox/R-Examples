plot.ITP2 <-
function(x,xrange=c(0,1),alpha1=0.05,alpha2=0.01,
                      ylab='Functional Data',main=NULL,lwd=1,col=c(1,2),pch=16,ylim=range(object$data.eval),
                      ...){
  object <- x
  if(alpha1 < alpha2){
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
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
    colors <- numeric(n)
    colors[which(object$labels==1)] <- col[1]
    colors[which(object$labels==2)] <- col[2]
    
    matplot(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=colors,lwd=lwd,ylim=ylim,...)
    
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
    colors <- numeric(n)
    colors[which(object$labels==1)] <- col[1]
    colors[which(object$labels==2)] <- col[2]
    
    matplot(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=colors,lwd=lwd,ylim=ylim,...)
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
    matplot(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=colors,lwd=lwd,add=TRUE,...)
    
    
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
    
  }else if(object$basis=='paFourier'){
    p <- length(object$pval_phase)
    J <- dim(object$data.eval)[2]
    n <- dim(object$data.eval)[1]
    xmin <- xrange[1]
    xmax <- xrange[2]
    abscissa.pval = 1:p
    Abscissa = seq(xmin,xmax,len=J)
    main.data <- paste(main,': Functional Data')
    main.data <- sub("^ : +", "", main.data)
    colors <- numeric(n)
    colors[which(object$labels==1)] <- col[1]
    colors[which(object$labels==2)] <- col[2]
    
    matplot(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=colors,lwd=lwd,ylim=ylim,...)
    
    ################################################################
    # pval phase
    main.p <- paste(main,': Adjusted p-values - phase')
    main.p <- sub("^ : +", "", main.p)
    plot(abscissa.pval,object$corrected.pval_phase,pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',xlab='Frequency')
    difference1 <- which(object$corrected.pval_phase<alpha1)
    if(length(difference1)>0){
      for(j in 1:length(difference1)){
        min.rect <- abscissa.pval[difference1[j]] - 0.5
        max.rect <- min.rect + 1
        rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = 'gray90',density=-2,border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = NULL,border='black')
    }
    difference2 <- which(object$corrected.pval_phase<alpha2)
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
    points(1:p,object$corrected.pval_phase,pch=pch)
    
    ################################################################
    # pval amplitude
    main.p <- paste(main,': Adjusted p-values - amplitude')
    main.p <- sub("^ : +", "", main.p)
    plot(abscissa.pval,object$corrected.pval_amplitude,pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',xlab='Frequency')
    difference1 <- which(object$corrected.pval_amplitude<alpha1)
    if(length(difference1)>0){
      for(j in 1:length(difference1)){
        min.rect <- abscissa.pval[difference1[j]] - 0.5
        max.rect <- min.rect + 1
        rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = 'gray90',density=-2,border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = NULL,border='black')
    }
    difference2 <- which(object$corrected.pval_amplitude<alpha2)
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
    points(1:p,object$corrected.pval_amplitude,pch=pch)
    
    
  }
  par(ask=FALSE)
}
