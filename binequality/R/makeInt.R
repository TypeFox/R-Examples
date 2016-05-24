makeInt <-
function(ints,hb, trans=NULL){
  up<-rep(ints[ ,2],times=hb)
  do<-rep(ints[,1], times=hb)
  do[which(do==0)]<-1.0001
  if(length(trans)>0){
    up<-trans(up)
    do<-trans(do)
  }#end if length(trans)
  
  out<-Surv(do,up,type='interval2')
  return(out)
}
