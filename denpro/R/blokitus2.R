blokitus2<-function(obj,blokki){
#
sar<-length(obj[1,]) #sarakkeiden maara 
riv<-length(obj[,1]) #rivien maara 
#
uusobj<-matrix(0,riv,sar+blokki)
uusobj[,1:sar]<-obj
#
return(uusobj)
}
