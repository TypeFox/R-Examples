legendg<-function(x,y=NULL,legend,fill=NULL,col=par("col"),          
 border=list("black"),lty,lwd,pch=NULL,angle=45,density=NULL,       
 bty="o",bg=par("bg"),box.lwd=par("lwd"),box.lty=par("lty"), 
 box.col=par("fg"),pt.bg=NA,cex=1,pt.cex=cex,pt.lwd=lwd,  
 pt.space=1,xjust=0,yjust=1,x.intersp=1,y.intersp=1,
 adj=c(0,0.5),text.width=NULL,text.col=par("col"),merge=FALSE,
 trace=FALSE,plot=TRUE,ncol=1,horiz=FALSE,title=NULL,
 inset=0,xpd,title.col=text.col) {

 if(missing(legend) && !is.null(y)) {
  legend<-y
  y<-NULL
 }
 if(is.list(x)) {
  y<-x$y
  x<-x$x
 }
 if(is.null(y)) {
  if(is.character(x)) {
   tablepos<-get.tablepos(x)
   x<-tablepos$x
   y<-tablepos$y
   xjust<-tablepos$xjust
   yjust<-tablepos$yjust
  }
 }
 if(!missing(xpd)) {
  oldxpd<-par("xpd")
  par(xpd=xpd)
 }
 legend.info<-legend(x=x,y=y,legend=legend,col=par("bg"),lty=1,       
  bty=bty,bg=bg,box.lwd=box.lwd,box.lty=box.lty, 
  box.col=par("fg"),pt.bg=NA,cex=cex,pt.cex=pt.cex,pt.lwd=pt.lwd,  
  xjust=xjust,yjust=yjust,x.intersp=x.intersp,y.intersp=y.intersp,
  adj=adj,text.width=text.width,text.col=text.col,merge=merge,
  trace=trace,plot=plot,ncol=ncol,horiz=horiz,title=title,
  inset=inset,title.col=title.col)
 if(!is.null(fill)) {
  if(length(border) < length(fill)) border<-rep(border,length(fill))
  # display the rectangles
  rectheight<-strheight("Q",cex=cex)
  if(length(adj) > 1) yadj<-adj[2] else yadj<-0.5
  for(nel in 1:length(fill)) {
   nrect<-length(fill[[nel]])
   rectspace<-(legend.info$text$x[nel]-legend.info$rect$left)
   lefts<-cumsum(c(legend.info$rect$left+rectspace*0.1,
    rep(0.8*rectspace/nrect,nrect-1)))
   rights<-lefts+0.7*rectspace/nrect
   bottoms<-rep(legend.info$text$y[nel]-yadj*rectheight,nrect)
   rect(lefts,bottoms,rights,bottoms+rectheight,col=fill[[nel]],
    border=ifelse(is.na(fill[[nel]]),NA,border[[nel]]))
  }
 }
 if(!is.null(pch)) {
  if(!is.list(col)) {
   mycol<-pch
   if(length(col) < length(mycol[[1]])) col<-rep(col,length.out=length(mycol[[1]]))
   for(nel in 1:length(col))
    mycol[[nel]]<-rep(col,length.out=length(mycol[[nel]]))
  }
  else mycol<-col
  lenpch<-length(pch)
  xright<-
   legend.info$text$x[1]-(legend.info$text$x[1]-legend.info$rect$left)*0.15
  for(nel in 1:lenpch) {
   npch<-length(pch[[nel]])
   pchwidth<-strwidth(pch[[nel]])*pt.space
   xpos<-rev((xright+pchwidth/2)-cumsum(pchwidth))
   ypos<-rep(legend.info$text$y[nel],npch)
   points(xpos,ypos,pch=pch[[nel]],col=mycol[[nel]],cex=pt.cex)
  }
 }
 if(!missing(xpd)) par(xpd=oldxpd)
 invisible(legend.info)
}
