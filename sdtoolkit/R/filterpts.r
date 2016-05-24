
filterpts <- function(dset, box, keepin=FALSE){

#box is assumed to be of "old" form, not Duong's.

#function will return dset, but with points removed as appropriate
#keepin=FALSE returns the points NOT covered by a box
#keepin=TRUE returns the points remaining inside a box.

  ndimr <- length(box[[1]])
  lmat <- matrix(nrow=nrow(dset),ncol=ndimr)

  for (i in 1:ndimr){
  
    lmat[,i] <- dset[,box[[1]][i]] > box[[2]][i,1] & dset[,box[[1]][i]] < box[[2]][i,2]
    
  }
  

  if(keepin){
  
    return(dset[apply(lmat,1,all),])
  
  }else{
  
    return(dset[!apply(lmat,1,all),])
    
  }
  
}
  
    