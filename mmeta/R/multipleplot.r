######################################################################################
######################################################################################
### Purpose: the wrapper of S3 method "plot" for object "multipletables"
### Input:   the object "multipletables"
### Output:  no return value.The plot may be saved as a pdf file or print in x11 windows
### Note:    the implement are functions "sideplot_multiple","overlapplot_multiple","forestplot"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
######################################################################################
plot.multipletables <- function(x,type=NULL,select=NULL,file=NULL, xlim=NULL,ylim=NULL,
                                xlabel=NULL,mar=NULL,xlog=TRUE,
                                addline=NULL,xlab=NULL,ylab=NULL,ciShow=TRUE,...) {
  if (!inherits(x, "multipletables"))
    stop("Use only with 'multiple' xs.\n")
  if (is.null(type)) stop("type is missing")
  measure<-x$measure
  if(xlog==TRUE & any(xlabel<0))
    stop("For log scale, label should be non-negative")                                              

  if (is.null(ylab)) ylab <- ""
  if(is.null(select)) {
    select <- seq(1:length(x$sample))
    if (type=="forest") select <- seq(1:(length(x$sample)+1))
  }
  if (type=="forest") if(max(select)>(length(x$sample)+1))
    stop("select is out of range. \n")
  if (type!="forest") if(max(select)>length(x$sample))
    stop("select is out of range. \n")
  if(type!="sidebyside" & type!="overlap" & type!="forest")
    stop("only 3 kinds of plots are available:sidebyside/overlap/forest")

  if(type=="sidebyside")
    sideplot_multiple(x,select=select,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,
                      xlabel=xlabel,addline=addline,file=file,mar=mar,...)
  if(type=="overlap")
    overlapplot_multiple(x,file=file,
                         xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,xlabel=xlabel,
                         addline=addline,select=select, mar=mar,...) 

  if(type=="forest") {
    forestplot(x,select=select,xlab=xlab,ylab=ylab,
               xlabel=xlabel,xlim=xlim,file=file,xlog=xlog,mar=mar,
               addline=addline,ciShow=ciShow,...)
  }
}

#########################################################################################
### Purpose: plot the forest plot
### Input:   S3 object "multipletables" 
### Output:  no return value.The plot may be saved as a pdf file or print in x11 windows
### Note:    the wrapper is function "plot.multipletables"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
#########################################################################################
forestplot <- function(object,select=NULL,xlab=NULL,ylab=NULL,
                       xlabel=xlabel,xlim=NULL,file=NULL,xlog=TRUE,mar=mar,
                       addline=NULL,ciShow=ciShow,...) {
  ##Quanlity control 
  measure<-object$measure
  report<-study_specifc(object)
      
  if (measure=="RD") {xlog<-FALSE; print("xlog is unavailable for RD")}
  j <- 1
  if(is.null(select)) select <- seq(1:(nrow(report[[j]])))
  nselect <- length(select)
  report.select <- report[[j]][select,]
  studynames <- rownames(report[[j]])[select]
  if (is.null(xlab)) {
    if (measure=="OR")  xlab <- "Odds ratio"
    if (measure=="RR")  xlab <- "Relative risk"
    if (measure=="RD")  xlab <- "Risk difference"
  }
  if (is.null(ylab)) ylab <- "Density"
  if (!is.null(xlim)) {
    xmax <- xlim[2]
    xmin <- xlim[1]
    truc.left <- which(report.select[,2]<xmin)
    truc.right <- which(report.select[,3]>xmax)
  }
  if (is.null(xlim)) {
    xmax <- max(report.select[,3])
    xmin <- min(report.select[,2])
  }
  if (is.null(xlabel)) {
    grid <- abs(xmax-xmin)/5
    mylabs <- round(c(xmin,xmin+grid,xmin+2*grid, xmin+3*grid,xmin+4*grid,xmax),2)
  }  else mylabs <- xlabel
  if (ciShow) {
    cishow <- vector()
    for (i in 1:nselect)
      cishow[i] <- paste(sprintf('%.1f',round(report.select[i,1],1))," (",
                         sprintf('%.1f',round(report.select[i,2],digits=1)),
                         ", ",sprintf('%.1f',round(report.select[i,3],1)),")",sep="")
  }
  if (xlog==TRUE) {
    report.select <- log10(report.select)
    xmax <- log10(xmax)
    xmin <- log10(xmin)
    if (!is.null(addline)) addline <- log10(addline)
    if (is.null(xlabel)) {
      grid <- abs(xmax-xmin)/5
      mylabs <- round(c(xmin,xmin+grid,xmin+2*grid, xmin+3*grid,xmin+4*grid,xmax),2)
    }  else mylabs <- log10(xlabel)
  }
  
  if (!is.null(file)) {
    origen.path <- getwd()
    savepath <- file.path(getwd(),"mmeta")
    dir.create(savepath,showWarnings = FALSE)
    setwd(savepath)
    pdf(paste(file,".pdf",sep=""))
    setwd(origen.path)
  }
 # if(is.null(file)) { dev.new() }
  if(!is.null(mar)) par(mar=mar) else par(mar=c(4, 5, 3, 6))
  plot(0,0,type="n", xlab=xlab, ylab=ylab, yaxt="n", xaxt="n", xaxs="i",
       yaxs="r",ylim=c(0,nselect+1), xlim=c(xmin,xmax),...)
  segments(report.select[nselect:1,2], 1:nselect ,report.select[nselect:1,3],1:nselect,lwd=2)
  segments(report.select[nselect:1,1], 1:nselect - 0.05, report.select[nselect:1,1], 1:nselect + 0.05, lwd=2)
  segments(report.select[nselect:1,2], 1:nselect - 0.05, report.select[nselect:1,2], 1:nselect+ 0.05, lwd=2)
  segments(report.select[nselect:1,3], 1:nselect - 0.05, report.select[nselect:1,3], 1:nselect + 0.05, lwd=2)
  u <- par("usr")
  axis(side=2, at=nselect:1, labels=studynames[1:nselect],las=1,cex.axis=0.8)
  if (!is.null(xlim))
    if (length(truc.left)>0)
      arrows(report.select[truc.left,1],nselect-truc.left+1,
             x1=rep(xmin,length(truc.left)),y1=nselect-truc.left+1,length = 0.15, angle = 15)
  if (!is.null(xlim))
    if(length(truc.right)>0)
      arrows(report.select[truc.right,1],nselect-truc.right+1,
             x1=rep(xmax,length(truc.right)),y1=nselect-truc.right+1,length = 0.15, angle = 15)
  
  if(ciShow) axis(side=4, at=nselect:1, labels=cishow[1:nselect],las=1,cex.axis=0.8)
  if(xlog==TRUE) axis(side=1, at=mylabs,labels=round(10^mylabs,1) )
  if(xlog==FALSE) axis(side=1, at=mylabs,labels=round(mylabs,1) )
  
  if (!is.null(addline)) {
    abline(v=addline,lty=2,col="blue")
    if(xlog==TRUE)
      if(!(addline%in%mylabs)) axis(side=1, at=addline,labels=10^addline)
    if(xlog==FALSE)
      if(!(addline%in%mylabs)) axis(side=1, at=addline,labels=addline)
  }

  if (!is.null(file)) {
    dev.off()
    cat( file, ".pdf have been saved in:", savepath, fill=TRUE)
  }
}
 
#############################################################################################
### Purpose: plot the posterir distribution in a side by side manner
### Input:   S3 object "multipletables" ,other input refers to the help file "plot.multiple"
### Output:  no return value.The plot may be saved as a pdf file or print in x11 windows
### Note:    the wrapper is function "plot.multipletables"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
#############################################################################################  
sideplot_multiple <- function(object,select=NULL,
                              xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,file=NULL,
                              xlabel=NULL,mar=mar,addline=NULL,...) {
  measure <- object$measure
  if(is.null(select)) select <- seq(1:length(object$sample))
  alpha <- object$alpha
  studynames.select <- object$studyname[select]
  density.select <- sample.select <- list()
  n.select <- length(select)
  sample.select <- object$sample[select]
  density.select <- object$density[select]
  priordens <- object$priordens

  if (is.null(xlab)) {
    if (measure=="OR")  xlab <- "Odds ratio"
    if (measure=="RR")  xlab <- "Relative risk"
    if (measure=="RD")  xlab <- "Risk difference"
  }
  if (is.null(ylab)) ylab <- "Density"
  k <- ceiling(n.select/4)
  if (n.select<3) rownumber<-1 else rownumber <- 2
  if (n.select>1) colnumber<- 2  else colnumber <- 1
  for (j in 1:k) {
    if (!is.null(file)) {
      origen.path <- getwd()
      savepath <- file.path(getwd(),"mmeta")
      dir.create(savepath,showWarnings = FALSE)
      filename <- paste(file,j,".pdf",sep="")
      setwd(savepath)
      pdf(filename)
      setwd(origen.path)
    }  else  { dev.new() }
    par(mfrow=c(rownumber,colnumber))
    if(!is.null(mar)) par(mar=mar) else par(mar=c(4, 3, 2, 1))

    ibegin <- (j-1)*4+1; iend<- min(j*4,n.select)
    for (i in ibegin:iend) {
      xmin <- quantile(sample.select[[i]],probs=alpha/10 ,na.rm=TRUE)
      xmax <- quantile(sample.select[[i]],probs=1-alpha/10 ,na.rm=TRUE)
      ymax <- max(density.select[[i]]$y[!is.na(density.select[[i]]$y)])*1.2
      ymin <- 0
      if (!is.null(xlim)) {
        xmax <- xlim[2]; xmin <- xlim[1]
      }
      if (!is.null(ylim)) {
        ymax <- ylim[2]; ymin <- ylim[1]
      }
      if(!is.null(xlabel)) {
        plot(density.select[[i]]$x, density.select[[i]]$y,
             type="l", lwd=2, xaxt="n", ylab=ylab, xlab=xlab,
             ylim=c(ymin,ymax),xlim=c(xmin,xmax),lty=1,...)
        axis(side=1, at=xlabel,labels=xlabel)
      }	else plot(density.select[[i]]$x, density.select[[i]]$y,
                  type="l", lwd=2,ylab=ylab, xlab=xlab,
                  ylim=c(ymin,ymax),xlim=c(xmin,xmax),lty=1,...)

      points(priordens$x,priordens$y,type='l', lty=3, lwd=2,...)

      if (!is.null(addline)) {
        abline(v=addline,lty=2,col="blue")
        if(!(addline%in%xlabel)) axis(side=1, at=addline,labels=addline )
      }
      legend("topright", c(studynames.select[i],"Prior"), lty=c(1,3), lwd=2, bty="n")
    }
    if (!is.null(file)){
      dev.off()
      cat( filename,"have been saved in:", savepath, fill=TRUE)
    }
  }
}      


#############################################################################################
### Purpose: plot the overlaid posterir distributions
### Input:   S3 object "multipletables" ,other input refers to the help file "plot.multiple"
### Output:  no return value.The plot may be saved as a pdf file or print in x11 windows
### Note:    the wrapper is function "plot.multipletables"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
############################################################################################## 
overlapplot_multiple <- function(object,file=NULL,
                                 xlab=NULL,ylab=NULL,select=NULL,xlim=NULL,ylim=NULL,
                                 xlabel=xlabel, addline=NULL, mar=mar,...) {
  measure <- object$measure
  if(is.null(select)) select <- seq(1:length(object$sample))
  
  alpha <- object$alpha
  studynames.select <- object$studyname[select]
  density.select <- sample.select <- list()
  n.select <- length(select)
  if(n.select<2) cat("select should be greater than 2","\n")

  sample.select <- object$sample[select]
  density.select <- object$density[select]

  if (is.null(xlab)) {
    if (measure=="OR")  xlab <- "Odds ratio"
    if (measure=="RR")  xlab <- "Relative risk"
    if (measure=="RD")  xlab <- "Risk difference"
  }
  if (is.null(ylab)) ylab <- "Density" 

  ##determine default range
  xmin <- ymax <- xmax <- vector()
  for(i in 1:n.select) {
    xmin[i] <- quantile(sample.select[[i]],probs=0.0001,na.rm=TRUE)
    xmax[i] <- quantile(sample.select[[i]],probs=0.9999 ,na.rm=TRUE)
    ymax[i] <- max(density.select[[i]]$y[!is.na(density.select[[i]]$y)])
  }
  if (is.null(xlim)) {
    xmin <- max(xmin)
    xmax <- min(xmax)
  } else {
    xmin<-xlim[1]; xmax<-xlim[2]
  }
  if (is.null(ylim)) {
    ymax <- max(ymax)
    ymin <- 0
  } else {
    ymin<-ylim[1]; ymax<-ylim[2]
  }

  if (!is.null(file))  {
    origen.path <- getwd()
    savepath <- file.path(getwd(),"mmeta")
    dir.create(savepath,showWarnings = FALSE)
    setwd(savepath)
    pdf(paste(file,".pdf",sep=""))
    setwd(origen.path)
  }
 # if(is.null(file)) { dev.new() }
   if(!is.null(mar)) par(mar=mar) else par(mar=c(4, 3, 3, 3))
  if(!is.null(xlabel)){
    plot(density.select[[1]]$x,density.select[[1]]$y, type="l",
         lwd=2,lty=1,xaxt="n",ylab=ylab, xlab=xlab,
         ylim=c(ymin,ymax),xlim=c(xmin,xmax),...)
    axis(side=1, at=xlabel,labels=xlabel)
  } else plot(density.select[[1]]$x, density.select[[1]]$y,
              type="l", lwd=2,lty=1,ylab=ylab, xlab=xlab,
              ylim=c(ymin,ymax),xlim=c(xmin,xmax),...)

  if (!is.null(addline)) {
    abline(v=addline,lty=2,col="blue")
    if(!(addline%in%xlabel)) axis(side=1, at=addline,labels=addline )
  }

  for(i in 2:n.select)
    points(density.select[[i]]$x, density.select[[i]]$y, type='l', lty=i, lwd=2,...)
  legend("topright", studynames.select, lty=1:n.select, lwd=rep(3,4),bty = "n")

  if (!is.null(file)) {
    dev.off()
    cat( file,".pdf have been saved in:", savepath,fill=TRUE)
  }
}

                                                                                                                                                                                     
                                                                                                                        

