crosshair <- function(x, ...) UseMethod("crosshair")

crosshair.meta4diag = function(x, est.type="mean", add=FALSE, save=F, xlim, ylim, xlab, ylab, ...){
  if(class(x)!="meta4diag"){stop("Wrong input given!")}
  if(!(est.type %in% c("mean","median","mode"))){stop("Argument \"est.type\" should be \"mean\",\"median\" or \"mode\".")}

  fitname = x$names.fitted
  fullname = paste("summary.fitted.(",fitname,")",sep="")
  
  if(est.type=="mean"){
    est.A = x[[fullname[1]]][,1]
    est.B = x[[fullname[2]]][,1]
  }
  if(est.type=="median"){
    est.A = x[[fullname[1]]][,4]
    est.B = x[[fullname[2]]][,4]
  }
  
  
  lb.A = x[[fullname[1]]][,3]
  ub.A = x[[fullname[1]]][,5]
  lb.B = x[[fullname[2]]][,3]
  ub.B = x[[fullname[2]]][,5]
  
  if(missing(xlim)){
    if(x$misc$model.type==1){xlim = c(1,0)}
    if(x$misc$model.type==2){xlim = c(0,1)}
    if(x$misc$model.type==3){xlim = c(1,0)}
    if(x$misc$model.type==4){xlim = c(0,1)}
  }else{
    if(x$misc$model.type==1){if(xlim[1]<xlim[2]){xlim = c(xlim[2],xlim[1])}}
    if(x$misc$model.type==2){if(xlim[2]<xlim[1]){xlim = c(xlim[2],xlim[1])}}
    if(x$misc$model.type==3){if(xlim[1]<xlim[2]){xlim = c(xlim[2],xlim[1])}}
    if(x$misc$model.type==4){if(xlim[2]<xlim[1]){xlim = c(xlim[2],xlim[1])}}
  }
  if(missing(ylim)){
    if(x$misc$model.type==1){ylim = c(0,1)}
    if(x$misc$model.type==2){ylim = c(0,1)}
    if(x$misc$model.type==3){ylim = c(1,0)}
    if(x$misc$model.type==4){ylim = c(1,0)}
  }else{
    if(x$misc$model.type==1){if(ylim[2]<ylim[1]){ylim = c(ylim[2],ylim[1])}}
    if(x$misc$model.type==2){if(ylim[2]<ylim[1]){ylim = c(ylim[2],ylim[1])}}
    if(x$misc$model.type==3){if(ylim[1]<ylim[2]){ylim = c(ylim[2],ylim[1])}}
    if(x$misc$model.type==4){if(ylim[1]<ylim[2]){ylim = c(ylim[2],ylim[1])}}
  }
  
  if(add){
    points(est.B, est.A, ...)
    arrows(lb.B, est.A, ub.B, est.A, angle=90, code=3, length=0.05, ...)
    arrows(est.B, lb.A, est.B, ub.A, angle=90, code=3, length=0.05, ...)
  }else{
    if(missing(xlab)){xlab=fitname[2]}
    if(missing(ylab)){ylab=fitname[1]}
    par(mfrow=c(1,1),mar=c(5.1, 4.1, 4.1, 2.1))
    plot(est.B, est.A, xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab, type="p",
         xaxs = "r",family="sans",xaxt="s",yaxt="s",bty="o",...)
    arrows(lb.B, est.A, ub.B, est.A, angle=90, code=3, length=0.05, ...)
    arrows(est.B, lb.A, est.B, ub.A, angle=90, code=3, length=0.05, ...)
  }
}