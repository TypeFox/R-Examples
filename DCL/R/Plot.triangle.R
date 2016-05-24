Plot.triangle<-function(triangle,Histogram=FALSE,tit='')
{
  Xtriangle<-as.matrix(triangle)
  m<-nrow(Xtriangle)
  
  # traditional plots
  if (Histogram==FALSE)
  {
 
  par(mfrow=c(2,1))
  par(mar=c(3,1.5,1.5,1.5),oma=c(2,0.5,0.5,0.2),mgp=c(1.5,0.5,0))
  # c(bottom, left, top, right)
  par(cex.main=1,cex.axis=0.9)
  
  matplot(t(Xtriangle),xlab='development period',ylab='',type='b',
          main=paste(tit, '-incremental'))
  # now plot the cummulative triangles
  Xtriangle.0<-Xtriangle
  Xtriangle.0[is.na(Xtriangle)]<-0
  cumX<-t(apply(Xtriangle.0,1,cumsum))
  cumX[is.na(Xtriangle)]<-NA
  matplot(t(cumX),xlab='development period',ylab='',type='b',
          main=paste(tit, '-cumulative'))
  par(mfrow=c(1,1))
  } else{
 
  colnames(Xtriangle)<-as.character(1:m)
  rownames(Xtriangle)<-as.character(1:m)
  cloud((Xtriangle), panel.3d.cloud = panel.3dbars,
        main=tit,xbase = 0.9, ybase = 0.9, 
        zlim = c(0, max(Xtriangle,na.rm=T)),
        scales = list(arrows = FALSE, just = "right"), 
        xlab = 'i', ylab ='j',zlab='',
        col.facet = level.colors(Xtriangle, 
                                 at = do.breaks(range(Xtriangle,na.rm=T), 20),
                                 col.regions = terrain.colors,
                                 colors = TRUE),
        colorkey = list(col = terrain.colors, 
                        at = do.breaks(range(Xtriangle,na.rm=T), 20)),
        screen = list(z = 40, x = -40))
  }
}
