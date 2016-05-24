blokitus.delt<-function(obj,blokki){
#
if (dim(t(obj))[1]==1) k<-1 else k<-length(obj[,1]) #rivien maara 
if (k==1){
  len<-length(obj)
  uusobj<-matrix(0,len+blokki,1)
  uusobj[1:len]<-obj
}
else{
  lev<-length(obj[1,])
  uusobj<-matrix(0,k+blokki,lev)
  uusobj[1:k,]<-obj
}
return(uusobj)
}
