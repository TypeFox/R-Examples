#-----------internal function ------------------
pvx.matrix<-function(theta_v,thres,xm_v=NULL){
  # func. by joerg-henrik heine jhheine(at)googlemail.com
  # ein dimension dazu und zum merken 
  # theta_v: ein vector oder zahl; 
  # thres: thurstonian thresholds eines items
  # xm_v: vector welche kategorie prob jeweils ausgegeben werden soll
  # korrigierte formel aus markus buch seite 330
  s<-0:length(thres)
  thres0<-c(0,thres)
  oben_v <- exp(apply((s%o%theta_v),2,function(x){x-cumsum(thres0)})) # ok - für theta_v als vector oder zahl
  unten_v<- apply( exp(apply((s%o%theta_v),2,function(x){x-cumsum(thres0)})) , 2 ,sum) # ok - für theta_v als vector oder zahl  
  px_v <- mapply(FUN=function(o,u){  o / u }, o=as.list(as.data.frame(oben_v)), u=as.list(unten_v) ) # u as list etc --> korrigiert am 10-3-2015 # ok - für theta_v als vector oder zahl 
  rownames(px_v)<-paste("cat",0:(length(thres0)-1),sep=".") 
  colnames(px_v)<-theta_v 
  P_v <- apply(px_v,2,sum) # ok - für theta_v als vector oder zahl
  if(length(xm_v)==0){return( (px_v) )}
  if(length(xm_v)!=0){mapply(function(p,ic){p[ic]}, as.list(as.data.frame(px_v)), xm_v)}
}     
#---------------------------------------------------  

