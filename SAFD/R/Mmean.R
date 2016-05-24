Mmean <-
function(XX,pic=0){
  #calculates the Minkowski mean of k=length(XX) polygonal fuzzy numbers with same levels
  #if necessary just use translator first to assure same alpha levels
  #if list contain missing values if necessary just use omitNA first
  k<-length(XX)
  X1<-Msum(XX)
  if(is.null(X1)==0){
   X2<-sc_mult(X1,1/k)
  
   #start possible plotting---------------------------------------------------
   if(pic==1){
      #calculate plot limits:
       lower<-rep(0,k)
       upper<-lower
        for (j in 1:k){
         lower[j]<-min(XX[[j]]$x)
         upper[j]<-max(XX[[j]]$x)
        }
       limx<-c(min(lower)-0.25,max(upper)+0.25)
     plot(XX[[1]],type="l", xlim=limx,lwd=1,xlab=NA, ylab=expression(alpha),cex.main=1, 
          main=paste("Sample and sample mean (in red)",sep=""))
     for (j in 2:k){
      lines(XX[[j]],type="l",lwd=1)
      }
     lines(X2,type="l",lwd=3,col="red")

    }
   #end possible plotting------------------------------------------------------
  invisible(X2)
  }
}
