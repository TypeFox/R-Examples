"qua.regressCOP.draw" <-
function(f=seq(0.1, 0.9, by=0.1),
         fs=0.5, cop=NULL, para=NULL, ploton=TRUE,
         wrtV=FALSE, col=c(4,2), lwd=c(1,2), lty=1, ...) {
  if(ploton) {
    plot(c(0,1),c(0,1), type="n", lwd=3,
         xlab="U, NONEXCEEDANCE PROBABILITY", ylab="V, NONEXCEEDANCE PROBABILITY")
  }
  Fs <- sort(c(f,fs))
  for(af in Fs) {
    mycol <- col[1]; mylwd <- lwd[1]
    if(af == fs) { mycol <- col[2]; mylwd <- lwd[2] }
    if(wrtV == FALSE) {
      R <- qua.regressCOP(f=af, cop=cop, para=para, ...)
    } else {
      R <- qua.regressCOP2(f=af, cop=cop, para=para, ...)
    }
    lines(R$U,R$V, col=mycol, lty=lty, lwd=mylwd)
  }
}
