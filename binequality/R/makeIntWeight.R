makeIntWeight <-
function(ints,hb, trans=NULL){
  up<-ints[ ,2]
  do<-ints[,1]
  #write.csv(up,file='up.csv',row.names=FALSE)
  #write.csv(do, file='do.csv',row.names=FALSE)
  do[which(do==0)]<-1.0001
  if(length(trans)>0){
    up<-trans(up)
    do<-trans(do)
  }#end if length(trans)
  
  int<-Surv(do,up,type='interval2')
  w<-hb/sum(hb)
  out<-list("int"=int,"weights"=w)
  return(out)
}
