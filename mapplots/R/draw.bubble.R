draw.bubble <-
function(x,y,z,maxradius=1,...){
  cex <- 2*maxradius/par('cxy')[2]/.375
  maxz <- max(z,na.rm=T)
  points(x,y,cex=cex*sqrt(z)/sqrt(maxz),...)
}

