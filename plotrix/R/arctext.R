arctext<-function(x,center=c(0,0),radius=1,start=NA,middle=pi/2,end=NA,stretch=1,
 clockwise=TRUE,cex=NA, ...) {

 oldcex <- par("cex")
 # have to do this to get strwidth to work
 if(is.na(cex)) cex <- oldcex
 par(cex=cex)
 xvec <- strsplit(x, "")[[1]]
 lenx <- length(xvec)
 xwidths <- stretch * strwidth(xvec)
 charangles <- xwidths/radius
 changrang <- range(charangles)
 # Make really narrow characters wider
 charangles[charangles < changrang[2]/2] <- changrang[2]/2
 if (!is.na(end)) start <- end + sum(charangles)
 if (is.na(start)) {
  if(clockwise) start <- middle + sum(charangles)/2
  else start <- middle - sum(charangles)/2
 }
 if(clockwise) {
  charstart <- c(start, start - cumsum(charangles)[-lenx])
  charpos <- charstart - charangles/2
 }
 else {
  charstart <- c(start, start + cumsum(charangles)[-lenx])
  charpos <- charstart + charangles/2
 }
 xylim <- par("usr")
 plotdim <- par("pin")
 ymult <- (xylim[4] - xylim[3])/(xylim[2] - xylim[1]) * plotdim[1]/plotdim[2]
 for (xchar in 1:lenx) {
  srt <- 180 * charpos[xchar]/pi - 90
  text(center[1] + radius * cos(charpos[xchar]), 
   center[2] + radius * sin(charpos[xchar]) * ymult,xvec[xchar], 
    adj=c(0.5, 0.5),srt=srt + 180 * !clockwise, ...)
 }
 par(cex=oldcex)
}
