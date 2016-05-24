"level.curvesCOP2" <-
function(cop=NULL, para=NULL, ploton=TRUE, lines=TRUE,
         plotMW=FALSE, ramp=TRUE, delv=0.001, delt=0.10,
         getlevel=NULL, ...) {

  if(! is.null(delt)) {
     if(delt <= 0 | delt > 0.5) {
        warning("Invalid 'delt' argument, must be (0,0.5] or NULL, setting to NULL")
        delt <- NULL
     }
  }
  if(! is.null(getlevel)) {
     if(getlevel < 0 | getlevel > 1) {
        warning("Invalid 'getlevel' argument, must be [0,1] or NULL, setting to NULL")
        getlevel <- NULL
     }
  }

  z <- list(level=getlevel, U=NULL, V=NULL)
  if(is.null(delt)) {
     Ts <- sort(getlevel)
  } else {
     Ts <- sort(c(getlevel, seq(0+delt,1-delt,delt)))
  }
  if(is.null(Ts)) {
     warning("The level set is NULL, check 'delt' and 'getlevel' arguments, returning NA")
     return(NA)
  }

  if(ploton) {
    plot(c(0,1), c(0,1), type="n",
         xlab="U, NONEXCEEDANCE PROBABILITY", ylab="V, NONEXCEEDANCE PROBABILITY")
  }

  for(t in Ts) {  # for each level t
    v <- c(t, t+delv/5, t+delv/2, seq(t+delv, 1-delv, by=delv), 1-delv/2, 1-delv/5, 1)
    v <- sort(v)
    u <- sapply(1:length(v), function(i) { COPinv2(cop=cop, v[i], t, para=para) } )
    if(lines) {
      if(ramp) {
        lines(u,v, lwd=(0.5+2*t), ...)
      } else {
        lines(u,v, ...)
      }
    }
    if(! is.null(getlevel) &
       isTRUE(all.equal(getlevel,t)) )
               z <- list(level=t, U=u, V=v)
  }
  if(plotMW) {
    abline(0,  1, lty=2, cex=0.5)
    abline(1, -1, lty=2, cex=0.5)
  }
  if(! is.null(getlevel)) return(z)
}
