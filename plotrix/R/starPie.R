starPie<-function(x,y,radext,values,maxval=NA,border=par("fg"),col=NA,
 prop.area=FALSE,label="",labelpos=1,radlab=NA,radprop=1.05,radlabcex=1) {

 valdim<-dim(values)
 if(is.null(valdim)) {
  npies<-1
  values<-matrix(values,nrow=1)
  valdim<-dim(values)
 }
 nfaces<-valdim[2]
 if(is.na(x[1])) x<-rep(1:valdim[1],valdim[2])
 if(is.na(y[1])) y<-rep(1:valdim[2],each=valdim[1])
 if(length(labelpos) < valdim[1])
  labelpos<-rep(labelpos,length.out=valdim[1]*valdim[2])
 if(is.na(col[1])) col<-rainbow(nfaces)
 if(length(col)<nfaces) col<-rep(col,length.out=nfaces)
 if(is.na(maxval[1])) maxval<-max(values)
 # get the y adjustment
 ymult<-getYmult()
 if(prop.area) values<-sqrt(values)
 for(spie in 1:valdim[1]) {
  angles<-5*pi/2-seq(0,pi*2,length.out=nfaces+1)
  labangles<-angles+diff(angles[1:2]/2)
  # for proportional area sectors, use the square root 
  # adjust the maximum value to the radius
  facerad<-radext*values[spie,]/maxval
  for(face in 1:nfaces) {
   xpos<-c(x[spie]+cos(angles[face])*facerad[face],x[spie],
    x[spie]+cos(angles[face+1])*facerad[face])
   ypos<-c(y[spie]+sin(angles[face])*facerad[face]*ymult,y[spie],
    y[spie]+sin(angles[face+1])*facerad[face]*ymult)
   polygon(xpos,ypos,col=col[face])
   segments(x[spie],y[spie],x[spie]+cos(angles[face])*radext,
    y[spie]+sin(angles[face])*radext*ymult)
   if(!is.na(radlab[1])) {
    par(xpd=TRUE)
    text(x[spie]+cos(labangles[face])*radext*radprop,
     y[spie]+sin(labangles[face])*radext*radprop,
     radlab[face],cex=radlabcex)
    par(xpd=FALSE)
   }
  }
  if(nchar(label[spie])) {
   x[spie]<-ifelse(labelpos[spie]%%2,x[spie],x[spie]+radext*(labelpos[spie]-3))
   y[spie]<-ifelse(labelpos[spie]%%2,y[spie]+radext*(labelpos[spie]-2),y[spie])
   hadj<-ifelse(labelpos[spie]%%2,0.5,labelpos[spie]==2)
   vadj<-ifelse(labelpos[spie]%%2,labelpos[spie]==3,0.5)
   text(x[spie],y[spie],label[spie],adj=c(hadj,vadj))
  }
 }
}
