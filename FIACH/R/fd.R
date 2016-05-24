fd<-function(input){
  if(is.character(input)==TRUE){
  rp<-read.table(input)}else{
    rp<-input
  }
  disp<-rp[,1:3]
  rad<-rp[,4:6]
  angDisp<-sin(rad)*50
  
  totDisp<-cbind(disp,angDisp)
  dif<-apply(totDisp,2,diff)
  dif.sq<-dif^2
  fd<-sqrt(rowSums(dif.sq))

  return(c(0,fd))
}

