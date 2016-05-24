"wrap.tor" <-
function(data,wrapav=TRUE,avestruc=NULL){
  
  wrap180 <- function(x) {
    x[which(x > 180, arr.ind=TRUE)] <-
      x[which(x > 180, arr.ind=TRUE)] - 360
    x[which(x < -180, arr.ind=TRUE)] <-
      x[which(x < -180, arr.ind=TRUE)] + 360
    x
  }
  
  if(!wrapav && is.null(avestruc))
    stop("Average structure is missing")
  if(is.vector(data)) {
    data <- matrix(data,ncol=1)
    return.vec = TRUE
  } else {
    return.vec = FALSE
  }

  
  avestruc.i<-avestruc
  datawrap <- NULL
  for(i in 1:ncol(data)) {
    struc <- data[,i]
    if(all(is.na(struc))) { struc <- rep(NA, length(struc)) } else {
    
      if(wrapav){ avestruc <- wrap180( mean(struc, na.rm=TRUE) ) } else {
        avestruc<-avestruc.i[i]
      }
	           
      difvar <- avestruc - struc
      while (length(difvar[ as.vector(na.omit(abs(difvar)>180)) ]) > 0) {
          
        struc[which(difvar > 180, arr.ind=TRUE)] <-
          struc[which(difvar > 180, arr.ind=TRUE)] + 360
          
        struc[which(difvar < -180, arr.ind=TRUE)] <-
          struc[which(difvar < -180, arr.ind=TRUE)] - 360
	    
        if(wrapav){avestruc <- wrap180( mean(struc, na.rm=TRUE) %% 360 )}
        
        difvar <- avestruc - struc
      }

    }
    datawrap <- cbind(datawrap,struc)
  }
  if(return.vec) datawrap = as.vector(datawrap)

  return(datawrap)
}
