plot.logistf.profile <-
function(x, type="profile", max1=TRUE, colmain="black", colimp="gray", plotmain=T, ylim=NULL, ...){
  obj<-x
 if(any(attr(obj,"class")[1]=="logistf.CLIP.profile")){
   x<-obj$beta
   if(type=="profile") {
     y<-obj$profile
     y.imp<-obj$profile.matrix
     y.label<-"Relative log profile penalized likelihood"
     } else
     if(type=="cdf") {
      y<-obj$cdf
      y.imp<-obj$cdf.matrix
      y.label<-"Cumulative distribution function"
     }   else
      if(type=="density") {
       x<-obj$beta[-length(obj$beta)]+diff(obj$beta)/2
       xy<-density(x, weights=diff(obj$cdf)/sum(diff(obj$cdf)))
       x<-xy$x
       y<-xy$y
#       y<-y/(sum(y)*sum(diff(obj$beta)))
# norm to AUC=1
       auc<-function(y,x){
        y.mean<-(y[-length(y)]+y[-1])/2
        x.diff<-diff(x)
        y.mean %*% x.diff
       }
       if(max1) y<-y/max(y)
       else y<-y/auc(y,x)
       ab<-t(apply(obj$cdf.matrix,1,diff))
       ab[ab<0]<-0
       ab<-ab/rowSums(ab)
       y.imp<-unlist(lapply(1:nrow(obj$cdf.matrix),function(Z) density(obj$beta[-length(obj$beta)]+diff(obj$beta)/2,
        weights=ab[Z,])$y))
       imputations<-nrow(obj$cdf.matrix)
       y.imp<-t(matrix(y.imp,length(y.imp)/imputations,imputations))
       if(max1) y.imp<-t(apply(y.imp,1,function(Z) Z/max(Z)))
       else y.imp<-t(apply(y.imp,1,function(Z) Z/auc(Z,x)))
       y.label<-"Posterior density"
      }
   plot(x=x, y=y, xlab=expression(gamma), ylab=y.label, type="l", lwd=2, col=colmain, ylim=ylim)
   if(type=="profile") abline(-3.841,0, lty=3)
   else if(type=="cdf") {
     abline(0.025,0, lty=3)
     abline(0.975,0, lty=3)
     }
   if(!is.null(y.imp)) {
     for(i in 1:nrow(y.imp)) lines(x=x, y=y.imp[i,], lwd=0.5, lty=2, col=colimp)
     if(plotmain) lines(x=x, y=y, lty=1, lwd=2, col=colmain)
     }
   } else
   if(any(attr(obj,"class")=="logistf.profile")) {
    x<-obj$beta
    if(type=="profile") {
      y<-obj$profile
      y.label<-"Relative log profile penalized likelihood"
      } else
      if(type=="cdf") {
       y<-obj$cdf
       y.label<-"Cumulative distribution function"
      } else
      if(type=="density") {
        x<-obj$beta[-length(obj$beta)]+diff(obj$beta)/2
        xy<-density(x, weights=diff(obj$cdf)/sum(diff(obj$cdf)))
        x<-xy$x
        y<-xy$y
#       y<-y/(sum(y)*sum(diff(obj$beta)))
# norm to AUC=1
        auc<-function(y,x){
         y.mean<-(y[-length(y)]+y[-1])/2
         x.diff<-diff(x)
         y.mean %*% x.diff
        }
        y<-y/auc(y,x)
        y.label<-"Posterior density"
       }
    plot(x=x, y=y, xlab=expression(gamma), ylab=y.label, type="l", lwd=2, col=colmain)
    if(type=="profile") abline(max(y)-3.841,0, lty=3)
    else if (type=="cdf") {
      abline(0.025,0, lty=3)
      abline(0.975,0, lty=3)
      }

   } else stop("Object must be of class profile or pool.profile.\n")
  


}

