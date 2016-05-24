"level.curvesCOP" <-
function(cop=NULL, para=NULL, ploton=TRUE, lines=TRUE,
         plotMW=FALSE, ramp=TRUE, delu=0.001, delt=0.10,
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
     tmp <- NULL # see comment below for seq() error trapping, it seems unlikely
     # that the trapping is needed here, but let us mimic the behavior.
     try(tmp <- seq(0+delt, 1-delt, by=delt), silent=FALSE)
     if(is.null(tmp)) try(tmp <- seq(0+delt, 1-delt, by=-delt), silent=FALSE)
     Ts <- sort(c(getlevel, tmp))
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
    # Concerning the tmp handling, there seems to be cases in which empirical copulas
    # are being used with small sample sizes that can trigger a 'wrong sign in 'by' argument
    # error on the following sequence, if so, let us try reversing the sequence and if that
    # fails then bail out enterly.
    tmp <- NULL
    try(tmp <- seq(t+delu, 1-delu, by=delu), silent=FALSE)
    if(is.null(tmp)) {
      warning("trying to compensate for 'by' error by reversing")
      try(tmp <- seq(t+delu, 1-delu, by=-delu), silent=FALSE)
    }
    u <- c(t, t+delu/5, t+delu/2, tmp, 1-delu/2, 1-delu/5, 1)
    u <- sort(u)
    v <- sapply(1:length(u), function(i) { COPinv(cop=cop, u[i], t, para=para) })
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
