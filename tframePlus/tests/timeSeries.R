  require("tframePlus")
  require("tfplot")
  if(!require("timeSeries")) stop("timeSeries not available, tests failed.")    


 Sys.info()

  #x11()
#  postscript(file="lite.out.ps",  paper="letter", horizontal=F, onefile=T)
#             # width=6, height=8, pointsize=10,


  all.ok <-  TRUE
  cat("tframe timeSeries test 1 ... ")
  z <- timeSeries(rnorm(100), 
     as.POSIXct(Sys.time() + sort(round(runif(100)*1e8)), "GMT")) 
  seriesNames(z) <- "random example"
  ok <- all(seriesNames(z) == c("random example"))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("unvariate seriesNames failed!\n") 
  ok <- is.tframed(z) & (inherits(z, "timeSeries"))
  all.ok <- ok
  if (ok) cat("ok\n") else cat("failed!\n") 


  cat("tframe timeSeries test 2 ... ")
  plot(z)
  tfplot(z)
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("failed!\n") 

  cat("tframe timeSeries test 3 ... ")
  y <- rnorm(100)
  tframe(y) <- tframe(z)  
  ok <- is.tframed(y) & (inherits(y, "timeSeries")) & (start(y) == start(z))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("failed!\n") 

  cat("tframe timeSeries test 4 ... ")
  ok <- all(tframe(y) == tframe(z))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("failed!\n") 


  cat("tframe timeSeries test 5 ... ")
  #  irregular series at random observation times
  # original z and the generated rnorm() series here do not have observations
  # at the same time, so the result (probably) has 220 Tobs.
  z <- tbind(z, timeSeries(rnorm(120),
          as.POSIXct(Sys.time() + sort(round(runif(120)*1e8)), "GMT")))
  seriesNames(z) <- c("random 1", "random 2")
  ok <- all(seriesNames(z) == c("random 1", "random 2"))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("multivariate seriesNames failed!\n") 
  # Tobs(z) will be 220 except in the random case of two equal time stamps
  y <- rnorm(Tobs(z))  
  tframe(y) <- tframe(z)  
  ok <- all(tframe(y) == tframe(z))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("failed!\n") 

  cat("tframe timeSeries test 6... ")
  plot(z)
  tfplot(z)
  ok <- is.tframed(z) & (inherits(z, "timeSeries"))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("failed!\n") 

  cat("tframe timeSeries test 7 ... ")
  y <- rnorm(Tobs(z)) 
  tframe(y) <- tframe(z)  
  ok <- is.tframed(y) & (inherits(y, "timeSeries")) & (start(y) == start(z))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("failed!\n") 

  cat("tframe timeSeries test 8 ... ")
  ok <- all(tframe(y) == tframe(z))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("failed!\n") 

  cat("tframe timeSeries test 9... ")
  tfplot(z, start=start(z))
  cat("ok\n") 


 cat("All tframe timeSeries tests completed")
 if (all.ok) cat(" OK\n") else cat(", some FAILED!\n") 
 
