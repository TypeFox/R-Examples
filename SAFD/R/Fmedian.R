Fmedian <-
function(XX,nl=101,pic=1){
  #calculates the 1-norm median of k=length(XX) polygonal fuzzy numbers with same levels
  #it is necessary use translator function first to work with a sufficiently large number of alpha levels
  m<-Mmean(XX)
  k<-length(XX)
  if(is.null(m)==0){
   X<-matrix(0,nrow=2*nl,ncol=k)
  YY<-vector("list",length=k)
 for(i in 1:k){
	YY[[i]]<-translator(XX[[i]],nl=nl)
      X[,i]<-YY[[i]]$x
   }
    Median<-c()
  for(j in 1:(2*nl)){
     Median[j]<-as.numeric(median(X[j,1:k]))
   }
   #return list containing the median
   Fmedian<-vector("list",length=1)
 levels<-seq(0,1,length=nl)
   Fmedian<-data.frame(x=Median,alpha=c(levels,rev(levels)))
   if(pic==1){
    #calculate plot limits:
    lower<-rep(0,k)
    upper<-lower
    for (i in 1:k){
      lower[i]<-min(XX[[i]]$x)
      upper[i]<-max(XX[[i]]$x)
     }
    limx<-c(min(lower)-0.25,max(upper)+0.25)
     plot(XX[[1]],type="l", xlim=limx,xlab=NA, ylab=expression(alpha),cex.main=1,lwd=1,
          main=paste("Sample, sample mean (blue) and sample median (red)",sep=""))
     for (i in 2:k){
      lines(XX[[i]],type="l",lwd=1)
      }
     lines(m,type="l", lwd=2.5,col="blue")
     lines(Fmedian,type="l", lwd=3,col="red")
    }
   #end possible plotting---------
   invisible(Fmedian)
  }
}
