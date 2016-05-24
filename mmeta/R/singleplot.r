##############################################################################################
### Purpose: the wrapper of S3 method "plot" for object "singletable"
### Input:   the object "singletable"
### Output:  no return value.The plot may be saved as a pdf file or be printed in x11 windows
### Note:    This function calls functions "sideplot_single", "overlapplot_single"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
##############################################################################################
plot.singletable <- function(x,type=type,file=NULL,select=c(1,2),
                             xlab=NULL,ylab=NULL,addline=NULL,xlim=NULL,
                             ylim=NULL,...) {
  if (!inherits(x, "singletable"))
    stop("Use only with 'singletable' objects.\n")

  if(max(select)>2) stop("select is out of range. \n")
  measure <- x$measure

  if(type!="sidebyside" & type!="overlap")
    stop("only 2 kinds of plots are available:sidebyside/overlap")
     
  if(type=="sidebyside")
    sideplot_single(x,file=file,
                    select=select,xlab=xlab,ylab=ylab,xlim=xlim,
                    ylim=ylim,addline=addline,...)
  if(type=="overlap")
    overlapplot_single(x,
                       file=file,select=select,xlab=xlab,ylab=ylab,xlim=xlim,
                       ylim=ylim,addline=addline,...)
}

################################################################################################
### Purpose: Plot the posterior distributuions in a side by side manner
### Input:   S3 object "singletable" ,other input refers to the help file "plot.single"
### Output:  no return value. The plot may be saved as a pdf file or be printed in x11 windows
### Note:    This function is called by wrapper function "plot.singletable"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
################################################################################################
sideplot_single <- function(object,file=NULL,
                            select=c(1,2),xlab=NULL,ylab=NULL,xlim=NULL,
                            ylim=NULL,addline=NULL,...) {
  alpha <- object$alpha
  measure <- object$measure
  nselect <- length(select)
  studynames.select <- object$studynames[select]
  j <- length(select)
  sample.select <- object$sample[select]
  density.select <- object$density[select]

  if (is.null(xlab)){
    if (measure=="OR")  xlab <- "Odds ratio"
    if (measure=="RR")  xlab <- "Relative risk"
    if (measure=="RD")  xlab <- "Risk difference"
  }

  if (is.null(ylab)) ylab <- "Density"

  if(!is.null(file)) {
    origen.path <- getwd()
    savepath <- file.path(getwd(),"mmeta")
    dir.create(savepath,showWarnings = FALSE)
    filename <- paste(file,".pdf",sep="")
    setwd(savepath)
    pdf(filename)
    setwd(origen.path)
  }
 # if(is.null(file)) { dev.new() }
  par(mfrow=c(1,j))

  for (i in 1:length(select)) {
    if (!is.null(xlim)) {
      xmax <- xlim[2]; xmin <- xlim[1]
    } else {
      xmin <- quantile(density.select[[i]]$x[density.select[[i]]$x!=0],probs=0.001 ,na.rm=TRUE)
      xmax <- quantile(density.select[[i]]$x[density.select[[i]]$x!=0],probs=0.999 ,na.rm=TRUE)
    }
    if (!is.null(ylim)) {
      ymax <- ylim[2]; ymin <- ylim[1]
    } else {
      ymin <- 0
      ymax <- max(density.select[[i]]$y[density.select[[i]]$y!=0])*1.2
    }
    plot(density.select[[i]]$x, density.select[[i]]$y, type="l", lwd=2,lty=1,
         axes = TRUE ,ylab=ylab, xlab=xlab,xlim=c(xmin,xmax),ylim=c(ymin,ymax),...)
    if (!is.null(addline)) {
      abline(v=addline,lty=3,col="blue")
      axis(side=1, at=addline,labels=addline)
    }
    legend("topright", studynames.select[i], lty=1, lwd=2,bty = "n")
  }
  if(!is.null(file)) {
    dev.off()
    cat(file," have been saved in:", savepath, fill=TRUE)
  }
}                                                                      
 
################################################################################################
### Purpose: Plot the overlaid posterior distributuions
### Input:   S3 object "singletable" ,other input refers to the help file "plot.single"
### Output:  no return value. The plot may be saved as a pdf file or be printed in x11 windows
### Note:    This function is called by function "plot.singletable"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
################################################################################################
overlapplot_single <- function(object,select=c(1,2),
                               file=NULL,xlab=NULL,ylab=NULL,xlim=NULL,
                               ylim=NULL,addline=NULL,...) {
  alpha <- object$alpha
  measure <- object$measure
  select <- c(1,2)
  nselect <- length(select)
  studynames.select <- object$studynames[select]
  method <- object$method
  sample.select <- object$sample[select]
  density.select <- object$density[select]
  gen.select <- object$gen[select]
  j <- length(select)

  if (is.null(xlab)) {
    if (measure=="OR") xlab <- "Odds ratio"
    if (measure=="RR") xlab <- "Relative risk"
    if (measure=="RD") xlab <- "Risk difference"
  }
  if (is.null(ylab)) ylab <- "Density"

  xmin <- xmax <- ymin <- ymax <- vector()
  for(i in 1:2) {
    xmin[i] <- quantile(density.select[[i]]$x[density.select[[i]]$x!=0],probs=0.001 ,na.rm=TRUE)
    xmax[i] <- quantile(density.select[[i]]$x[density.select[[i]]$x!=0],probs=0.999 ,na.rm=TRUE)
    ymin[i] <- 0
    ymax[i] <- max(density.select[[i]]$y)*1.2
  }

  if (is.null(xlim)) {
    xmax <- xmax[1]
    xmin <- xmin[1]
  } else {
    xmax <- xlim[2]
    xmin <- xlim[1]
  }

  if (is.null(ylim)) {
    ymax <- max(ymax)
    ymin <- 0
  } else {
    ymax <- ylim[2]
    ymin <- ylim[1]
  }

  if(!is.null(file)) {
    origen.path <- getwd()
    savepath <- file.path(getwd(),"mmeta")
    dir.create(savepath,showWarnings = FALSE)
    filename <- paste(file,".pdf",sep="")
    setwd(savepath)
    pdf(filename)
    setwd(origen.path)
  }
  #if(is.null(file)) { dev.new() }
  plot(density.select[[1]]$x, density.select[[1]]$y, type="l", lwd=2, axes = TRUE,
       ylim=c(ymin,ymax), xlim=c(xmin,xmax),ylab=ylab, xlab=xlab,...)
  for(i in 1:2) {
    points(density.select[[i]]$x, density.select[[i]]$y, type='l', lty=i, lwd=2)
    if (!is.null(addline)) {
      abline(v=addline,lty=3,col="blue")
      axis(side=1, at=addline,labels=addline )
    }
  }
  legend("topright",  studynames.select, lty=1:2, lwd=rep(3,4),bty = "n")
  if(!is.null(file)) {
    dev.off()
    cat(file, "have been saved in:", savepath,fill=TRUE)
  }
}
