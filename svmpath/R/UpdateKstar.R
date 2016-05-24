UpdateKstar<-function(Kstar,Kell,Krow,y){
  if(length(y)==1){
    advec<-c(y,Krow)
    Kstar<-rbind(Kstar,advec)
    cbind(Kstar,c(advec,Kell))
  }
  else{
    adrect<-cbind(y,Krow)
    Kstar<-rbind(Kstar,adrect)
    cbind(Kstar,rbind(t(adrect),Kell))
  }
}
  
