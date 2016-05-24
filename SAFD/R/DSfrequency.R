DSfrequency <-
function(XX,IV=c(0,1),pic=1,nl=101){
 #XX fuzzy sample (list)
 #IV is the interval for which the (levelwise) frequency will be calculated
 #interpretation: levelwise Dempster-Shafer frequencies
 #nl is the number of levels for which the levelwise Dempster-Shafer frequencies are calculated
 #by default nl=101
 k<-length(XX)
 temp<-rep(0,k)
  for (i in 1:k){
   temp[i]<-checking2(XX[[i]],0)
  }
 if(min(temp)==0){
   print("One or more elements of the input data defines no polygonal fuzzy number")
   }
 if(min(temp)==1){
  X<-translator(XX[[1]],nl)
  YY<-vector("list",length=k)
  YY[[1]]<-X
 
  hitting<-rep(0,nl)
  contained<-rep(0,nl)
  if(X$x[1]<=IV[2]&X$x[2*nl]>=IV[1]){
   hitting<-ifelse(X$x[1:nl]>IV[2]|X$x[(2*nl):(nl+1)]<IV[1],0,1)
   contained<-ifelse(X$x[1:nl]>=IV[1]&X$x[(2*nl):(nl+1)]<IV[2],1,0)
   }
  if(k==1){
   f<-data.frame(x=c(contained,hitting),alpha=X$alpha)
   if(pic==1){
     plot(YY[[1]],type="l", lwd=0.1,xlab=NA, ylab=expression(alpha),cex.main=1, col="gray50",
           main=paste("Sample and chosen interval",sep=""))
     Left<-data.frame(x=rep(IV[1],2),alpha=c(-0.05,1.05))
     Right<-data.frame(x=rep(IV[2],2),alpha=c(-0.05,1.05))
     lines(Left,type="l",col="red",lwd=2)
     lines(Right,type="l",col="red",lwd=2)

     dev.new()
     plot(f,type="l",xlim=c(0,1),ylim=c(0,1),
     main=paste("Levelwise Dempster Shafer frequency of the interval [",IV[1], ",",IV[2],"]",sep=""),
       cex.main=1,xlab=NA, ylab=expression(alpha))
       lines(f,type="p",cex=0.2)
     }
   invisible(f)
   }
  if(k>1){
   for (i in 2:k){
    X<-translator(XX[[i]],nl)
    YY[[i]]<-X
    hitting_dazu<-rep(0,nl)
    contained_dazu<-rep(0,nl)
    if(X$x[1]<=IV[2]&X$x[2*nl]>=IV[1]){
     hitting_dazu<-ifelse(X$x[1:nl]>IV[2]|X$x[(2*nl):(nl+1)]<IV[1],0,1)
     contained_dazu<-ifelse(X$x[1:nl]>=IV[1]&X$x[(2*nl):(nl+1)]<IV[2],1,0)
    }
    hitting<-hitting+hitting_dazu
    contained<-contained+contained_dazu
    }
   hitting<-1/k*hitting
   contained<-1/k*contained
  f<-data.frame(x=c(contained,hitting[nl:1]),alpha=YY[[1]]$alpha)
  if(pic==1){
      lower<-rep(0,k)
       upper<-lower
        for (j in 1:k){
         lower[j]<-min(YY[[j]]$x)
         upper[j]<-max(YY[[j]]$x)
        }
       limx<-c(min(lower),max(upper))
     plot(YY[[1]],type="l", xlim=limx,lwd=0.1,xlab=NA, ylab=expression(alpha),cex.main=1, col="gray50",
          main=paste("Sample and chosen interval",sep=""))
      for (j in 2:k){
      lines(YY[[j]],type="l",lwd=0.1,col="gray50")
      }
      Left<-data.frame(x=rep(IV[1],2),alpha=c(-0.05,1.05))
      Right<-data.frame(x=rep(IV[2],2),alpha=c(-0.05,1.05))
      lines(Left,type="l",col="red",lwd=2)
      lines(Right,type="l",col="red",lwd=2)
      dev.new()
      plot(f,type="l",xlim=c(0,1),ylim=c(0,1),
      main=paste("Levelwise Dempster Shafer frequency of the interval [",IV[1], ",",IV[2],"]",sep=""),
      cex.main=1,xlab=NA, ylab=expression(alpha))
      lines(f,type="p",cex=0.2)      }
  invisible(f)
 }
 }
}
