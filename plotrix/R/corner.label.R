corner.label<-function(label=NULL,x=-1,y=1,xoff=NA,yoff=NA,figcorner=FALSE,...) {

 if(is.na(xoff)) xoff<-strwidth("m")/2
 if(is.na(yoff)) yoff<-strheight("m")/2
 par.usr<-par("usr")
 xpos<-par.usr[(3+x)/2]
 ypos<-par.usr[(3+y)/2+2]
 if(figcorner) {
  par.pin<-par("pin")
  xplotrange<-par.usr[2]-par.usr[1]
  yplotrange<-par.usr[4]-par.usr[3]
  par.mai<-par("mai")
  xmar<-xplotrange*par.mai[3+x]/par.pin[1]
  ymar<-yplotrange*par.mai[2+y]/par.pin[2]
  xpos<-xpos+x*xmar
  ypos<-ypos+y*ymar
 }
 if(!is.null(label)) {
  if(figcorner) par(xpd=TRUE)
  text(xpos-x*xoff,ypos-y*yoff,label,adj=c((1+x)/2,(1+y)/2),...)
  if(figcorner) par(xpd=FALSE)
 }
 return(list(x=xpos,y=ypos))
}
