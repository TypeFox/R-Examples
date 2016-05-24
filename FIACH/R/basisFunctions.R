basisFunctions <-function(RT,orth=TRUE){
  can.hrf<- spmHrf(RT)
  dp<-1
  p<-can.hrf$p
  p[6]<-p[6]+dp
  D<- (can.hrf$hrf - spmHrf(RT,p)$hrf)/dp
  can.hrf$hrf<-cbind(can.hrf$hrf, D)
  p[6]<- p[6] - dp
  
  dp<-0.01;
  p[3]<-p[3] +dp
  D<- (can.hrf$hrf[,1] - spmHrf(RT,p)$hrf)/dp
  can.hrf$hrf<-cbind(can.hrf$hrf, D)
  p[3]<- p[3] -dp
  
  if(orth==TRUE){ret<-spmOrth(can.hrf$hrf)}else{
    ret<-can.hrf$hrf
  }
  return(ret)
}

