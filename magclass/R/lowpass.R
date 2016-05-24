lowpass <- function(x,i=1, fix=NULL) {
  
  if(!is.null(fix)) warning("Fixing start or end does might modify the total sum of values! Use fix=NULL to let the total sum unchanged!")
  
  if(i==0) return(x)
  
  if(is.magpie(x)) {
    for(k in 1:dim(x)[1]) {
      for(j in if(is.null(getNames(x))) 1 else getNames(x)) {
        x[k,,j] <- lowpass(as.vector(x[k,,j]),i=i,fix=fix)
      }
    }  
  } else {
    l <- length(x)
    for(j in 1:i) {
      y <- x
      x[2:(l-1)] <- (y[1:(l-2)] + 2*y[2:(l-1)] + y[3:l])/4  
      if(is.null(fix)){
        x[1] <- (3*y[1]+y[2])/4
        x[l] <- (3*y[l]+y[l-1])/4
      }
      else if (fix=="start") x[l] <- (3*y[l]+y[l-1])/4
      else if (fix=="end") x[1] <- (3*y[1]+y[2])/4
      else if (fix!="both") stop(paste("Option \"",fix,"\" is not available for the \"fix\" argunemt!",sep=""))    }
  }  
  return(x)
}