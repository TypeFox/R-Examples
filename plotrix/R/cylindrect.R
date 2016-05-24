cylindrect<-function(xleft,ybottom,xright,ytop,col,border=NA,gradient="x",
 nslices=50) {
 
 rgbval<-col2rgb(col)/255
 maxrgb<-max(rgbval)
 cols<-matrix(0,nrow=3,ncol=6)
 for(i in 1:3) {
  if(rgbval[i] == maxrgb) delta<-1-rgbval[i]
  else delta<-(1-0.2*(maxrgb-rgbval[i])/sum(maxrgb-rgbval))-rgbval[i]
  cols[i,]<-c(rgbval[i]+0.3*delta,rgbval[i]+0.6*delta,rgbval[i]+0.9*delta,
   rgbval[i]+0.6*delta,rgbval[i]+0.3*delta,rgbval[i])
 }
 gradient.rect(xleft,ybottom,xright,ytop,cols[1,],cols[2,],cols[3,],
  gradient=gradient,nslices=nslices,border=border)
 invisible(cols)
}
