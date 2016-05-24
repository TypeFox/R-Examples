print.inputlist <-
  function (x,...) 
  {
   if (!inherits(x, "inputlist") | !is.list(x) | (length(x) != 10)) 
      stop(paste(sQuote("x"), "is not a correct input parameters list"))

cat('###  Summary of the input parameters: ###
     \n Initial time t0 = ',x$t0,
    '\n Initial state x0 = ',x$x0, 
    '\n Final time Tfin (ms) = ', x$Tfin,
    '\n Timestep delta t (ms) = ',x$deltat, 
    '\n Threshold type = ',x$Sty, 
    '\n Number of spikes to be simulated = ', x$M, '\n') 

## Parameter quadflag controls whether a full reconstruction of the 
# FPT density function via numerical quadrature is required
if (x$quadflag == 1) cat ('\n Numerical Integration up to Tfin required \n \n') 
}