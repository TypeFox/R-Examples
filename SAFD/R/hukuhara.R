hukuhara <-
function(X,Y,pic=0){
 #it is checked if the hukuhara difference Y-X exists and if yes it is calculated and returned
 #first check validity of input data
 temp_sum<-Msum(list(X,Y))
 if(is.null(temp_sum)==0){
  nl<-nrow(X)/2
  dif<-data.frame(x=Y$x-X$x,alpha=X$alpha)
  #calculate for each alpha-level the interval-hukuhara diff and check if this
  #family of intervals is decreasing in alpha ->if yes then return the polygonal fuzzy number that
  #is the hukuhara diff
  if(pic==1){
    plot(X,type="l",xlim=c(min(c(X$x,Y$x)),max(c(X$x,Y$x))),xlab=NA, ylab=expression(alpha))
    lines(Y,type="l")
    }
    a<-checking(dif,0)
  if(a==0){
   print("Hukuhara difference Y-X does not exist")
   }
 if(a==1){
  if(pic==1){
   plot(X,type="l",xlim=c(min(c(X$x,Y$x)),max(c(X$x,Y$x))),xlab=NA, ylab=expression(alpha),
        lwd=2, main=paste("Fuzzy numbers and their Hukuhara difference (in red)",sep=""),cex.main=1)
   lines(Y,type="l",lwd=2)
   lines(dif,type="l",lwd=3,col="red")
   }
  invisible(dif)
  }
  }
 }
