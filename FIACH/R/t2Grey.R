t2Grey<-function(B0,relax=TRUE){
  if(relax){return(1.74*B0+7.77)}else{
    return(1/(1.74*B0+7.77)*1000)
  }
}

