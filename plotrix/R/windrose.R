bin.wind.records<-function(winddir,windspeed,ndir=8,radians=FALSE,
 speed.breaks=c(0,10,20,30)) {
 # the windagg matrix has wind speeds as rows and wind directions as columns
 windagg<-matrix(0,ncol=ndir,nrow=5)
 dir.breaks<-rep(0,ndir)
 if(radians) {
  angleinc<-2*pi/ndir
  dir.breaks[1]<-pi/ndir
 }
 else {
  angleinc<-360/ndir
  dir.breaks[1]<-180/ndir
 }
 for(i in 2:ndir) dir.breaks[i]<-dir.breaks[i-1]+angleinc
 nspeeds<-length(speed.breaks)+1
 for(i in 1:length(winddir)) {
  dir<-1
  while(winddir[i] > dir.breaks[dir] && dir < ndir) dir<-dir+1
  if(winddir[i] > dir.breaks[ndir]) dir<-1
  speed<-1
  while(windspeed[i] > speed.breaks[speed] && speed < nspeeds) speed<-speed+1
  windagg[speed,dir]<-windagg[speed,dir]+1
 }
 windagg<-100*windagg/sum(windagg)
 return(windagg)
}

oz.windrose.legend<-function(maxpct=20,scale.factor=30,
 speed.col=c("#dab286","#fe9a66","#ce6733","#986434"),
 speed.width=NA,legend.pos=NA) {
  
 wdnames<-c("E","NE","N","NW","W","SW","S","SE")
 if(is.na(speed.width[1])) speed.width<-maxpct*1:4/100
 if(is.na(legend.pos[1])) legend.pos<-maxpct*1.25
 draw.circle(-maxpct/2,legend.pos,maxpct/20)
 angles<-seq(0,7*pi/4,by=pi/4)
 for(i in 1:8) {
  x<-cos(angles[i])*c(maxpct/20,maxpct/16,maxpct/10)-maxpct/2
  y<-sin(angles[i])*c(maxpct/20,maxpct/16,maxpct/10)+legend.pos
  segments(x[1],y[1],x[2],y[2])
  text(x[3],y[3],wdnames[i],cex=0.7)
 }
 wsnames<-c("1-10","10-20","20-30","30+")
 draw.circle(-maxpct/30,legend.pos,maxpct/30)
 for(i in 1:length(speed.col)) {
  x<-c(i-1,i)*maxpct/4
  y<-c(legend.pos-speed.width[i],legend.pos+speed.width[i])
  polygon(c(x[1],x[1],x[2],x[2]),c(y[1],y[2],y[2],y[1]),col=speed.col[i])
  text((x[1]+x[2])/2,legend.pos-maxpct/15,wsnames[i],cex=0.7)
 }
 text(-maxpct/30,legend.pos+maxpct/10,"Calm")
 text(maxpct/2,legend.pos+maxpct/10,"km/h")
 text(maxpct/2,legend.pos-maxpct/8,paste("Scale factor = ",scale.factor,"%",sep=""))
}

oz.windrose<-function(windagg,maxpct=20,wrmar=c(4,5,6,5),
 scale.factor=30,speed.col=c("#dab286","#fe9a66","#ce6733","#986434"),
 speed.width=NA,show.legend=TRUE,legend.pos=NA,...) {
 
 if(is.na(speed.width[1])) speed.width<-maxpct*1:4/100
 if(is.na(legend.pos[1])) legend.pos<-maxpct*1.25
 oldmar<-par("mar")
 par(mar=wrmar,xpd=TRUE)
 plot(0,xlim=c(-maxpct,maxpct),ylim=c(-maxpct,maxpct),type="n",
  axes=FALSE,xlab="",ylab="",...)
 winddim<-dim(windagg)
 calm.radius<-sum(windagg[1,])/winddim[2]
 rad<-10
 while(rad<=maxpct) {
  draw.circle(0,0,rad+calm.radius,border="gray")
  boxed.labels(rad+calm.radius,maxpct/10,paste(rad,"%",sep=""),ypad=0.7,border=FALSE)
  rad<-rad+10
 }
 draw.circle(0,0,calm.radius,border="gray")
 angle.inc<--2*pi/winddim[2]
 angles<-seq(5*pi/2,pi/2+angle.inc,by=angle.inc)
 descpct<-order(colSums(windagg),decreasing=TRUE)
 for(i in descpct) {
  next.radius<-calm.radius
  for(j in 2:winddim[1]) {
   xinner<-cos(angles[i])*next.radius
   xouter<-cos(angles[i])*(windagg[j,i]+next.radius)
   yinner<-sin(angles[i])*next.radius
   youter<-sin(angles[i])*(windagg[j,i]+next.radius)
   next.radius<-sqrt(xouter*xouter+youter*youter)
   # find the four points for each rectangle
   xoffset<-cos(angles[i]-pi/2)*speed.width[j-1]
   yoffset<-sin(angles[i]-pi/2)*speed.width[j-1]
   polygon(c(xinner-xoffset,xinner+xoffset,xouter+xoffset,xouter-xoffset),
    c(yinner-yoffset,yinner+yoffset,youter+yoffset,youter-yoffset),
    col=speed.col[j-1])
  }
 }
 text(-maxpct,maxpct/4,paste("Calm ",round(sum(windagg[1,]),1),"%",sep=""),col="blue")
 if(show.legend)
  oz.windrose.legend(maxpct=maxpct,scale.factor=scale.factor,speed.col=speed.col,
   speed.width=speed.width,legend.pos=legend.pos)
 par(oldmar)
}
