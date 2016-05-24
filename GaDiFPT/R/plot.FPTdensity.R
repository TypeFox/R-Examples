plot.FPTdensity <-
  function (x, ...) 
  {
   if (!inherits(x, "FPTdensity") | !is.list(x) | (length(x) != 3)) 
      stop(paste(sQuote("x"), "is not a correct FPTdensity object"))

   # hazard rate function lambda = g0/(1-gg0)
   lambda  <- x$g0/(1-x$gg0)
   N1p1 <- length(lambda)
   
   RStflag <- get("RStudioflag")
   
   # produces plots
   
   if(!RStflag) dev.new()
   plot(x$time,x$g0,type="l",lwd=1,main="FPT density",xlab="time(ms)",ylab="")
   if(!RStflag) dev.new()
   plot(x$time,x$gg0,type="l",lwd=1,main="FPT distribution",xlab="time(ms)",ylab="")
   if(!RStflag) dev.new()
   plot(x$time,lambda,type="l",lwd=1,main="FPT hazard function",xlab="time(ms)",ylab="")

}