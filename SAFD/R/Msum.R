Msum <-
function(XX,pic=0){
  #calculates the Minkowski sum of k polygonal fuzzy numbers with same levels
  #if necessary just use translator first to assure same alpha levels
  #if check is not ok or if sample elements have different alpha levels NA is returned
  #if list contain missing values if necessary just use omitNA first
  k<-length(XX)
  temp<-rep(0,k)
  for (i in 1:k){
   temp[i]<-checking(XX[[i]],0)
  }
  ok<-1
  if(min(temp)==0){
   print(paste("One or more elements of the list don't define a polygonal fuzzy number"))
   print(paste("use checking function to specify the violations in the definition of a poligonal fuzzy number"))
   ok<-0
   }
  if(min(temp)==1){
   if(length(XX)==1){R<-XX[[1]]}
   if(length(XX)>=2){
    number<-rep(0,k)
    for (i in 1:k){
     number[i]<-nrow(XX[[i]])
     }
    if(max(number)!=min(number)){
     print("use translator function to assure that list elements are compatible (same alpha levels)")
     ok<-0
     }
    if(max(number)==min(number)){
     equal<-rep(0,k-1)
      for (i in 1:(k-1)){
      equal[i]<-max(abs(XX[[1]]$alpha-XX[[i]]$alpha))
      }
     if(max(equal)>0){
      print("use translator function to assure that list elements are compatible (same alpha levels)")
      ok<-0
     }
     }
    }
    if(length(XX)>=2&ok==1){
      R<-XX[[1]]
       for (i in 2:k){
        R$x<-R$x+XX[[i]]$x
       }
    }
   }
   if(ok==1){
   #start possible plotting---------------------------------------------------
   if(pic==1){
      #calculate plot limits:
       lower<-rep(0,k)
       upper<-lower
        for (j in 1:k){
         lower[j]<-min(XX[[j]]$x)
         upper[j]<-max(XX[[j]]$x)
        }
       limx<-c(min(c(lower,R$x))-0.25,max(c(upper,R$x))+0.25)
     plot(XX[[1]],type="l",lwd=2, xlim=limx,xlab=NA, ylab=expression(alpha),cex.main=1,
          main=paste("Sum (in red) of fuzzy numbers",sep=""))
     for (j in 2:k){
      lines(XX[[j]],type="l",lwd=2)
      }
     lines(R,type="l",lwd=3,col="red")

    }
   #end possible plotting------------------------------------------------------

invisible(R)
  }
}
