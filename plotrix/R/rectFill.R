rectFill<-function(x1,y1,x2,y2,fg=par("fg"),bg=par("bg"),xinc=NA,yinc=NA,
 pch=1,pch.cex=1,pch.col=par("fg"),...) {

 if(is.na(xinc[1])) xinc<-strwidth(pch,cex=pch.cex)
 if(is.na(yinc[1])) yinc<-strheight(pch,cex=pch.cex)
 nrect<-length(x1)
 if(length(fg) < nrect) fg<-rep(fg,nrect)
 if(length(bg) < nrect) bg<-rep(bg,nrect)
 if(length(xinc) < nrect) xinc<-rep(xinc,nrect)
 if(length(yinc) < nrect) yinc<-rep(yinc,nrect)
 if(length(pch) < nrect) pch<-rep(pch,nrect)
 if(length(pch.cex) < nrect) pch.cex<-rep(pch.cex,nrect)
 if(length(pch.col) < nrect) pch.col<-rep(pch.col,nrect)
 for(frect in 1:nrect) {
  rect(x1[frect],y1[frect],x2[frect],y2[frect],col=bg[frect],border=fg[frect])
  if(yinc[frect] > 0) {
   xpos<-seq(x1[frect]+xinc[frect]/2,x2[frect]+xinc[frect]/2,by=xinc[frect])
   lenxpos<-length(xpos)
   ypos<-seq(y1[frect]+yinc[frect]/2,y2[frect]+yinc[frect],by=yinc[frect])
   lenypos<-length(ypos)
   xpos<-rep(xpos,each=lenypos)
   ypos<-rep(ypos,lenxpos)
   clip(x1[frect],y1[frect],x2[frect],y2[frect])
   points(xpos,ypos,pch=pch[frect],cex=pch.cex[frect],col=pch.col[frect],...)
  }
 }
 do.call(clip,as.list(par("usr")))
}
