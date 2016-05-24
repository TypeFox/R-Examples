# staxlab produces staggered axis tick labels
# note that barplot() tends to mess things up by plotting an X axis 
# even when axes=FALSE

staxlab<-function(side=1,at,labels,pos=NA,nlines=2,top.line=0.5,
 line.spacing=0.8,srt=NA,ticklen=0.03,adj=1,...) {

 if(missing(labels)) labels<-at
 nlabels<-length(labels)
 if(missing(at)) at<-1:nlabels
 axislim<-par("usr")[3:4-2*side%%2]
 if(any(at < axislim[1]) || any(at > axislim[2]))
  warning("Some axis labels are off the plot")
 if(is.na(srt)) {
  linepos<-rep(top.line,nlines)
  if(nlines > 1) for(i in 2:nlines) linepos[i]<-linepos[i-1]+line.spacing
  linepos<-rep(linepos,ceiling(nlabels/nlines))[1:nlabels]
  axis(side=side,at=at,labels=rep("",nlabels),pos=pos)
  mtext(text=labels,side=side,line=linepos,at=at,...)
 }
 else {
  linewidth<-strheight("M")
  xylim<-par("usr")
  if(side == 1) {
   xpos<-at
   if(par("ylog")) ypos<-10^(xylim[3]-ticklen*(xylim[4]-xylim[3]))
   else ypos<-xylim[3]-ticklen*(xylim[4]-xylim[3])-top.line*linewidth
  }
  if(side == 3) {
   xpos<-at
   if(par("ylog")) ypos<-10^(xylim[4]+ticklen*(xylim[4]-xylim[3]))
   else ypos<-xylim[4]+ticklen*(xylim[4]-xylim[3])+top.line*linewidth
  }
  if(side == 2) {
   ypos<-at
   if(par("xlog")) xpos<-10^(xylim[1]-ticklen*(xylim[2]-xylim[1]))
   else xpos<-xylim[1]-ticklen*(xylim[2]-xylim[1])-top.line*linewidth
  }
  if(side == 4) {
   ypos<-at
   if(par("xlog")) xpos<-10^(xylim[2]+ticklen*(xylim[2]-xylim[1]))
   else xpos<-xylim[2]+ticklen*(xylim[2]-xylim[1])+top.line*linewidth
  }
  par(xpd=TRUE)
  text(xpos,ypos,labels,srt=srt,adj=adj,...)
  par(xpd=FALSE)
 }
}
