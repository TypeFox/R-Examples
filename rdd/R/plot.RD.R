#' Plot of the Regression Discontinuity
#' 
#' Plot the relationship between the running variable and the outcome
#' 
#' @method plot RD
#' @param x \code{rd} object, typically the result of \code{\link{RDestimate}}
#' @param gran the granularity of the plot. This specifies the number of points to 
#' either side of the cutpoint for which the estimate is calculated.
#' @param bins if the dependent variable is binary, include the number of bins 
#' within which to average
#' @param which identifies which of the available plots to display. For a sharp 
#' design, the only possibility is \code{1}, the plot of the running variable 
#' against the outcome variable. For a fuzzy design, an additional plot, \code{2},
#' may also be displayed, showing the relationship between the running variable and
#' the treatment variable. Both plots may be displayed with \code{which=c(1,2)}.
#' @param range the range of values of the running variable for which to plot. This
#' should be a vector of length two of the format \code{c(min,max)}. To plot
#' from the minimum to the maximum value, simply enter \code{c("min","max")}. The default 
#' is a window 20 times wider than the first listed bandwidth from the \code{rd} object, 
#' truncated by the min/max values of the running variable from the data.
#' @param ... unused
#' @details It is important to note that this function will only plot the discontinuity
#' using the bandwidth which is first in the vector of bandwidths passed to \code{RDestimate}
#' @include RDestimate.R
#' @importFrom stats predict lm
#' @importFrom graphics lines plot
#' @importFrom grDevices dev.flush
#' @export
#' @author Drew Dimmery <\email{drewd@@nyu.edu}>


plot.RD <- function(x,gran=400,bins=100,which=1,range,...) {
  frm<-FALSE
  if("frame" %in% names(x$call)) frm<-eval.parent(x$call$frame)
  if(!frm){
    x$call$frame<-TRUE
    x$call$verbose<-FALSE
    x<-eval.parent(x$call)
  }
  d<-as.data.frame(x$frame)

  if(length(x$na.action)>0)
    d<-d[-x$na.action,]
  
  if("kernel" %in% names(x$call)) 
    kern<-eval.parent(x$call$kernel)
  else 
    kern<-"triangular"
  
  if("cutpoint" %in% names(x$call)) 
    cut<-eval.parent(x$call$cutpoint)
  else
    cut<-0
  
  bw<-x$bw[1]
  
  if(missing(range)) {
     range<-c(cut-10*bw,cut+10*bw)
		 if(range[1]<min(d$X)) range[1]<-min(d$X)
		 if(range[2]>max(d$X)) range[2]<-max(d$X)
  }
     
  if(range[1]=="min")
    range[1]<-min(d$X)
  if(range[2]=="max")
    range[2]<-max(d$X)
  range<-as.double(range)
  
  rdplot<-function(d) {
    d.l<-data.frame(X=d$X[d$X<cut],Y=d$Y[d$X<cut])
    lval<-seq(range[1],cut,length.out=(gran%/%2))
    lest<-vector(length=(gran%/%2))
    llwr<-vector(length=(gran%/%2))
    lupr<-vector(length=(gran%/%2))
    for(i in 1:(gran%/%2)) {
      sub<-d.l$X>=(lval[i]-bw) & d.l$X<=(lval[i]+bw)
      w<-kernelwts(X=d.l$X[sub],center=lval[i],bw=bw,kernel=kern)
      ly<-d.l$Y[sub]
      lx<-d.l$X[sub]
      if(length(lx)<=2)
        pred<-rep(NA,3)
      else
        pred<-predict(lm(ly~lx,weights=w),interval="confidence",newdata=data.frame(lx=lval[i]))
      lest[i]<-pred[1]
      llwr[i]<-pred[2]
      lupr[i]<-pred[3]
    }
    
    d.r<-data.frame(X=d$X[d$X>=cut],Y=d$Y[d$X>=cut])
    rval<-seq(cut,range[2],length.out=(gran%/%2))
    rest<-vector(length=(gran%/%2))
    rlwr<-vector(length=(gran%/%2))
    rupr<-vector(length=(gran%/%2))
    for(i in 1:(gran%/%2)) {
      sub<-d.r$X>=(rval[i]-bw) & d.r$X<=(rval[i]+bw)
      w<-kernelwts(X=d.r$X[sub],center=rval[i],bw=bw,kernel=kern)
      ry<-d.r$Y[sub]
      rx<-d.r$X[sub]
      if(length(rx)<=2)
        pred<-rep(NA,3)
      else
        pred<-predict(lm(ry~rx,weights=w),interval="confidence",newdata=data.frame(rx=rval[i]))
      rest[i]<-pred[1]
      rlwr[i]<-pred[2]
      rupr[i]<-pred[3]
    }
    
    #plot to the left
    if(length(unique(d$Y))==2) {
      #DO THIS for when the outcome is dichotomous
      ep<-(max(d$X)-min(d$X))/(2*bins)
      nX<-seq(min(d$X)-ep,max(d$X)+ep,length=bins+1)
      nY<-rep(NA,length(nX))
      for(i in (1:(length(nX)-1))){
        if(sum(!is.na(d$Y[d$X>nX[i] & d$X<=nX[i+1]]))==0)
          next
        nY[i]<-sum(d$Y[d$X>nX[i] & d$X<=nX[i+1]],na.rm=TRUE)/sum(!is.na(d$Y[d$X>nX[i] & d$X<=nX[i+1]]))
      }
      sub<-nX>=range[1] & nX<=range[2]
      subl<-lval>=range[1] & lval<=range[2]
      subr<-rval>=range[1] & rval<=range[2]
      plot(nX,nY,
           type="p",pch=20,cex=.5,col="black",
           xlim=c(range[1],range[2]),
           ylim=c(min(c(llwr[subl],rlwr[subr]),na.rm=T),
                  max(c(lupr[subl],rupr[subr]),na.rm=T)),
           xlab=NA,
           ylab=NA,
           main=NA
      )
    } else {
      subl<-lval>=range[1] & lval<=range[2]
      subr<-rval>=range[1] & rval<=range[2]
      plot(d$X,d$Y,
           type="p",pch=20,cex=.5,col="black",
           xlim=c(range[1],range[2]),
           ylim=c(min(c(llwr[subl],rlwr[subr]),na.rm=T),
                  max(c(lupr[subl],rupr[subr]),na.rm=T)),
           xlab=NA,
           ylab=NA,
           main=NA
      )
    } 
    #plot to the left
    lines(lval,lest,
          lty=1,lwd=2,col="black",type="l"
    )
    
    lines(lval,llwr,
          lty=2,lwd=1,col="black",type="l"
    )
    lines(lval,lupr,
          lty=2,lwd=1,col="black",type="l"
    )
    
    #plot to the right
    lines(rval,rest,
          lty=1,lwd=2,col="black",type="l"
    )
    lines(rval,rlwr,
          lty=2,lwd=1,col="black",type="l"
    )
    lines(rval,rupr,
          lty=2,lwd=1,col="black",type="l"
    )
  }
  if(x$type=="sharp" | 1%in%which){
    rdplot(d)
    dev.flush()
  }
  if(x$type=="fuzzy" & 2%in%which){
    d$Y<-d$Z
    rdplot(d)
    dev.flush()
  }
}
