"diagCOP" <-
function(cop=NULL, para=NULL, secondary=FALSE,
         ploton=TRUE, lines=TRUE, delt=0.005, ...) {
  if(ploton) {
    plot(c(0,1), c(0,1), type="n",
         xlab="U=V, NONEXCEEDANCE PROBABILITY",
         ylab="H, NONEXCEEDANCE PROBABILITY")
  }
  T <- seq(0+delt,1-delt,delt)
  C <- sapply(T, function(x) {
          	        y <- x
          	        if(secondary) y <- 1 - y
          	        return( cop(x,y, para=para, ...)) } )
  if(lines & ! is.null(dev.list())) lines(T,C, ...)
  txt <- "primary"
  if(secondary) txt <- "secondary"
  return(list(t=T, diagcop=C, diagtype=txt))
}
