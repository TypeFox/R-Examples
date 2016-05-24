Gmedian <-
function(behavior,phaseX,v){
  
  abmedian<-tapply(behavior, phaseX,median,na.rm=T)
  omedian<- abmedian[names(abmedian)==v]
  abline(h=omedian,col="black",lwd=3)
}
