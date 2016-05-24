ablineclip<-function(a=NULL,b=NULL,h=NULL,v=NULL,reg=NULL, 
 coef=NULL,untf=FALSE,x1=NULL,x2=NULL,y1=NULL,y2=NULL,...) {

 if(!is.null(c(x1,x2,y1,y2))) {
  oldclip<-par("usr")
  # if any clipping perimeters are not supplied, use the existing plot edges
  if(is.null(x1)) x1<-oldclip[1]
  if(is.null(x2)) x2<-oldclip[2]
  if(is.null(y1)) y1<-oldclip[3]
  if(is.null(y2)) y2<-oldclip[4]
  clip(x1,x2,y1,y2)
  abline(h=oldclip[4]+1)
  clip(x1,x2,y1,y2)
 }
 abline(a=a,b=b,h=h,v=v,reg=reg,coef=coef,untf=untf,...)
 if(!is.null(c(x1,x2,y1,y2))) do.call("clip",as.list(oldclip))
}
