PicCalc <-
function(data){
  data[is.na(data)]<-0
  a<-rowSums(data==1)
  b<-rowSums(data==-1)
  return(mean(1-(a^2+b^2)/(a+b)^2,na.rm=TRUE))
}
