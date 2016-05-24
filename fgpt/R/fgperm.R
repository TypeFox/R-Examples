fgperm <-
function(xy,z=1:dim(xy)[1], scale, group=1, iter=999, ratio=1, FUN=fyshuffle, ..., add.obs=FALSE, as.matrix=FALSE){  
  
  runrand <- function(t){
    
    m <- array(NA,dim=c(dim(xy)[1],9))
    m[,1] <- z
    m[,2] <- group
    
    theta <- runif(1,0,pi)
    
    m[,3] <- (xy[,1] * cos(theta) + xy[,2] * sin(theta))/(scale*sqrt(ratio))
    m[,4] <- (xy[,2] * cos(theta) - xy[,1] * sin(theta))/(scale/sqrt(ratio))
    m[,5] <- runif(1,min(m[,3])-1,min(m[,3]))
    m[,6] <- runif(1,min(m[,4])-1,min(m[,4]))
    m[,7:8] <- ceiling(m[,3:4] - m[,5:6])
    cells <- unique(m[,c(2,7,8)])  
    m[,9] <- sapply(1:dim(m)[1],function(i,x,cells){which(cells[,1]==x[i,1] & cells[,2]==x[i,2] & cells[,3]==x[i,3])}, cells=cells, x=m[,c(2,7,8)])
    
    return(ave(m[,1],m[,9],FUN=FUN))
  } 
  
  if(as.matrix==TRUE){
    output <- sapply(1:iter,runrand)    
  }else{
    output <- lapply(1:iter,runrand)
  }
  if(add.obs==TRUE){
    output <- c(observed=list(z),output)
  }
  return(output)
}
