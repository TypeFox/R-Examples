"orderPed"<-function(ped, time_born=NULL){
  
  if(length(time_born)!=0){
    reorder<-order(time_born)
  }else{
    reorder<-order(kindepth(ped[,1],ped[,2],ped[,3]), decreasing=FALSE)
  }

  ped[reorder,]
}
