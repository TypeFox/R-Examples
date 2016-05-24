"sectionCOP" <-
function(f, cop=NULL,  para=NULL, wrtV=FALSE, dercop=FALSE, delt=0.005,
              ploton=TRUE, lines=TRUE,  xlab="NONEXCEEDANCE PROBABILITY", ...) {

  if(wrtV) {
     message("Triggering Horizontal Section logic: v = constant")
     txt <- "horizontal section"
  } else {
     message("Triggering Vertical Section logic: u = constant")
     txt <- "vertical section"
  }

  if(ploton) plot(c(0,1), c(0,1), type="n", xlab=xlab, ...)

  T <- seq(0+delt,1-delt,delt)
  C <- vector(mode="numeric", length=length(T))
  if(dercop) {
    C <- sapply(T, function(x) {
                     ifelse(wrtV, return( derCOP2(x,f, cop=cop, para=para) ),
                                  return(  derCOP(f,x, cop=cop, para=para) )) } )
  } else {
    C <- sapply(T, function(x) {
                      ifelse(wrtV, return( cop(x,f, para=para) ),
                                   return( cop(f,x, para=para) )) } )
  }
  if(lines & ! is.null(dev.list())) lines(T,C, ...)
  return(list(t=T, seccop=C, wrt=txt, fvalue=f, isderivative=dercop))
}

