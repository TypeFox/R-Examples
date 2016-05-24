dirSeg<-function(x,deg=TRUE){
   dx<-x[,3]-x[,1]
   dy<-x[,4]-x[,2]
   direct<-ifelse(dx>=0 & dy >0,atan(dx/dy),ifelse(dx>0 & dy<=0, pi/2+atan(abs(dy)/dx),ifelse(dx<=0 & dy<0, pi+atan(abs(dx)/abs(dy)),1.5*pi+atan(dy/abs(dx)))))
   if (deg) direct<-direct/pi*180
   direct
}