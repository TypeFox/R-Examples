lgdDisplayCircles <-
function(posLegCircles,pt,varSize,var,txtCexLeg,
                            legendCircles,
                            circleCol,lgdRnd){
  inset <- c(0.01, 0.01)
  rmax<-max(pt[,varSize],na.rm=TRUE)
  rhalf<-sqrt(((rmax*rmax*pi)/3)/pi)
  rmin<-min(pt[,varSize],na.rm=TRUE)
  rLeg<-c(rmax,rhalf,rmin)
  vmax<-max(pt[,var],na.rm=TRUE)
  vhalf<-vmax/3
  vmin<-min(pt[,var],na.rm=TRUE)
  rVal<-c(vmax,vhalf,vmin)
  usr <- par("usr")
  inset.x <- (usr[2] - usr[1]) * inset[1]
  inset.y <- (usr[4] - usr[3]) * inset[2]
  
  circ.width <- rLeg[1]*2
  space.width <- rLeg[1]*0.1
  text.width <- max(strwidth(round(rVal,lgdRnd),cex=txtCexLeg))
  total.width <- max(c((circ.width + space.width + text.width),strwidth(legendCircles,cex=txtCexLeg)))
  total.height <- circ.width + strheight("lp",cex=txtCexLeg)
  
  left <- switch(posLegCircles, 
                 bottomright = , topright = , 
                 right = usr[2] - total.width - inset.x, 
                 bottomleft = ,left = , 
                 topleft = usr[1] + inset.x, 
                 bottom = , top = , 
                 center = (usr[1] + usr[2] - total.width)/2)
  
  top <- switch(posLegCircles, 
                bottomright = , bottom = , 
                bottomleft = usr[3] +  total.height + inset.y + strheight(rLeg[1],cex=txtCexLeg)/2 , 
                topleft = ,top = , 
                topright = usr[4] - inset.y -  strheight(legendCircles,cex=txtCexLeg), 
                left = ,right = , 
                center = (usr[3] + usr[4] + total.height)/2)
  text(left, top  , labels = legendCircles, cex = txtCexLeg, adj=c(0, 0))
  top<-top-strheight("lp",cex=txtCexLeg)
  symbols(rep(left+rLeg[1],3),y=top-(2*rLeg[1])+rLeg,circles=rLeg,add=TRUE,inches=FALSE,bg=circleCol,lwd=0.5)
  for (i in 1:3){
    segments (left+rLeg[1],
              top-(2*rLeg[1])+2*rLeg[i],
              left+(rLeg[1]*2.05),
              top-(2*rLeg[1])+2*rLeg[i],
              lwd=0.5)
  }
  text(x=rep(left+(rLeg[1]*2.05),3),y=top-(2*rLeg[1])+2*rLeg,round(rVal,lgdRnd),cex=txtCexLeg,srt=0,adj=0)
}
