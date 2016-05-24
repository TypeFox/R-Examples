SCOC <- function(times, y=NULL, parms, Flux, ...) {

  if (is.null(y)){
    meanDepo <- mean(approx(Flux[,1],Flux[,2], xout=seq(1,365,by=1))$y)
    y <- meanDepo/parms
  } else if (length(y) != 1)
    stop ("length of state variable vector should be 1")

  if (length(parms) != 1)
    stop ("length of parameter vector should be 1")

  names(y) <- c("C")

  out <- vode(y, times, func = "scocder",
    parms = parms, dllname = "deSolve",
    initforc="scocforc",  forcings=Flux,
    initfunc = "scocpar", nout = 2,
    outnames = c("Mineralisation","Depo"),...)
  out
}

