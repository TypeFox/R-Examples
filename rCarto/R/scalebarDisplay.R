scalebarDisplay <-
function(scalebar,scalebarSize,scalebarText,txtCexScalebar){
  if(scalebar==TRUE){
    l<-locator(n=1)
    SpatialPolygonsRescale(layout.scale.bar(),offset=c(l$x,l$y),
                           scale=scalebarSize,fill=c("black"),plot.grid=F)
    text(l$x+scalebarSize/2,(l$y),paste(scalebarText,"\n\n",sep=""),cex=txtCexScalebar)
  }
}
