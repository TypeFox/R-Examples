  require("tframePlus")
  require("tfplot")
  if(!require("tis")) stop("tis not available, tests failed.")    
  if(!require("zoo")) stop("zoo not available, tests failed.")    


 Sys.info()

  #x11()
#  postscript(file="lite.out.ps",  paper="letter", horizontal=F, onefile=T)
#             # width=6, height=8, pointsize=10,


  all.ok <-  TRUE
  cat("tframe tis test 1 ... ")
  z <- tis(rnorm(100), start="2007-01-01", frequency=52)
  seriesNames(z) <- "random example"
  ok <- all(seriesNames(z) == c("random example"))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("unvariate seriesNames failed!\n") 
  ok <- is.tframed(z) & (inherits(z, "tis"))
  all.ok <- ok
  if (ok) cat("ok\n") else cat("failed!\n") 


  cat("tframe tis test 2 ... ")
  plot(z)
  tfplot(z)
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("failed!\n") 

  cat("tframe tis test 3 ... ")
  y <- rnorm(100)
  tframe(y) <- tframe(z)  
  ok <- is.tframed(y) & (inherits(y, "tis")) & (start(y) == start(z))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("failed!\n") 

  cat("tframe tis test 4 ... ")
  ok <- all(tframe(y) == tframe(z))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("failed!\n") 


  cat("tframe tis test 5 ... ")
  #  tbind
  z <- tbind(z, tis(rnorm(100), start="2007-01-01", frequency=52))
  seriesNames(z) <- c("random 1", "random 2")
  ok <- all(seriesNames(z) == c("random 1", "random 2"))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("multivariate seriesNames failed!\n") 
  # Tobs(z) 100
  y <- rnorm(Tobs(z)) 
  tframe(y) <- tframe(z)  
  ok <- all(tframe(y) == tframe(z))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("failed!\n") 

  cat("tframe tis test 6... ")
  plot(z)
  tfplot(z)
  ok <- is.tframed(z) & (inherits(z, "tis"))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("failed!\n") 

  cat("tframe tis test 7 ... ")
  y <- rnorm(Tobs(z))
  tframe(y) <- tframe(z)  
  ok <- is.tframed(y) & (inherits(y, "tis")) & (start(y) == start(z))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("failed!\n") 

  cat("tframe tis test 8 ... ")
  ok <- all(tframe(y) == tframe(z))
  all.ok <- all.ok & ok 
  if (ok) cat("ok\n") else cat("failed!\n") 


 cat("All tframe tis tests completed")
 if (all.ok) cat(" OK\n") else cat(", some FAILED!\n") 
 
