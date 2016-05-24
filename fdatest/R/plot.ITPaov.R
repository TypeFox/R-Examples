plot.ITPaov <-
function(x,xrange=c(0,1),alpha1=0.05,alpha2=0.01,plot.adjpval=FALSE,
                      ylim=range(x$data.eval),col=1,
                      ylab='Functional Data',main=NULL,lwd=1,pch=16,
                      ...){
  if(alpha1 < alpha2){
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
  
  object <- x
  p <- length(object$pval.F)
  J <- dim(object$data.eval)[2]
  n <- dim(object$data.eval)[1]
  xmin <- xrange[1]
  xmax <- xrange[2]
  abscissa.pval = seq(xmin,xmax,len=p)
  Abscissa = seq(xmin,xmax,len=J)
  par(ask=T) 
  main.f <- paste(main,': Functional Data and F-test')
  main.f <- sub("^ : +", "", main.f)
  
  
  matplot(Abscissa,t(object$data.eval),type='l',col=col,main=main.f,ylab=ylab,ylim=ylim,lwd=lwd,...)
  difference1 <- which(object$corrected.pval.F < alpha1)
  if (length(difference1) > 0) {
    for (j in 1:length(difference1)) {
      min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
      max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
      rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray90", density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
  }
  
  difference2 <- which(object$corrected.pval.F < alpha2)
  if (length(difference2) > 0) {
    for (j in 1:length(difference2)) {
      min.rect <- abscissa.pval[difference2[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
      max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
      rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray80", density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
  }
  matplot(Abscissa,t(object$data.eval),type='l',col=col,add=TRUE,lwd=lwd,...)
  
  formula <- object$call$formula
  mf <- model.frame(formula)
  nvar <- dim(object$corrected.pval.factors)[1]
  names_all <- colnames(object$design.matrix)
  interaz <- grep(':',names_all)
  
  for(var in 1:(dim(object$corrected.pval.factors)[1])){
    var.name = rownames(object$corrected.pval.factors)[var]
    main.t <- paste(main,': factor',var.name,sep=' ')
    main.t <- sub("^ : +", "", main.t)
    
    if(length(grep(':',var.name))>0){ # sto plottando interazione
      var12 <- strsplit(var.name,':')
      var1 <- var12[[1]][1]
      var2 <- var12[[1]][2]
      dummy.test1 <- grep(var1,names_all)
      dummy.test2 <- grep(var2,names_all)
      dummy.test <- intersect(dummy.test1,dummy.test2)
      colors <- object$design.matrix[,dummy.test]
      if(length(dim(colors))>1){
        colors <- (apply(colors,1,paste,collapse=''))
      }
      colors <- as.factor(colors)
    }else{ #sto plottando un fattore
      dummy.test <- grep(var.name,names_all)
      dummy.test <- setdiff(dummy.test,interaz)
      colors <- object$design.matrix[,dummy.test]
      if(length(dim(colors))>1){
        colors <- (apply(colors,1,paste,collapse=''))
      }
      colors <- as.factor(colors)
    }
    
    
    matplot(Abscissa,t(object$data.eval),type='l',col=colors,ylim=ylim,lwd=1,main=main.t,ylab=ylab,...)
    difference1 <- which(object$corrected.pval.factors[var,] < alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray90", density = -2, border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
    }
    difference2 <- which(object$corrected.pval.factors[var,] < alpha2)
    if (length(difference2) > 0) {
      for (j in 1:length(difference2)) {
        min.rect <- abscissa.pval[difference2[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray80", density = -2, border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
    }
    matlines(Abscissa,t(object$data.eval),type='l',col=colors,...)
    #lines(ascissa,coeff.teo[1,],lty=2,add=TRUE,type='l',col=1,lwd=2)
    abline(h=0,lty=2,col=1)
  }
  #########################################################
  #plot of adjusted p-values
  if(plot.adjpval==TRUE){
    main.p <- paste(main,': Adjusted p-values - F-test')
    main.p <- sub("^ : +", "", main.p)
    Abscissa <- abscissa.pval
    plot(Abscissa,object$corrected.pval.F,pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',...)
    difference1 <- which(object$corrected.pval.F<alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray90", density = -2, border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
    }
    difference2 <- which(object$corrected.pval.F<alpha2)
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
    points(Abscissa,object$corrected.pval.F,pch=pch)
    
    for(var in 1:(dim(object$corrected.pval.factors)[1])){
      var.name = rownames(object$corrected.pval.factors)[var]
      main.p <- paste(main,': Adjusted p-values - factor',var.name)
      main.p <- sub("^ : +", "", main.p)
      plot(Abscissa,object$corrected.pval.factors[var,],pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',...)
      difference1 <- which(object$corrected.pval.factors[var,]<alpha1)
      if (length(difference1) > 0) {
        for (j in 1:length(difference1)) {
          min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
          max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
          rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray90", density = -2, border = NA)
        }
        rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
      }
      difference2 <- which(object$corrected.pval.factors[var,]<alpha2)
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
      points(Abscissa,object$corrected.pval.factors[var,],pch=pch)
    }
  }
  par(ask=F) 
}
