# A short set of tests of the tframe class methods. 
  require("tframe")

 Sys.info()

  #x11()
#  postscript(file="lite.out.ps",  paper="letter", horizontal=F, onefile=T)
#             # width=6, height=8, pointsize=10,


verbose <- TRUE
synopsis <- TRUE

  all.ok <-  TRUE
  if (synopsis & !verbose) cat("All tframe tests ...")
  if (verbose) cat("tframe test 1 ... ")
  #tspvector <- tframed(1:100, list(start=c(1981,3), frequency=4))
  tspvector <- ts(1:100, start=c(1981,3), frequency=4)
  dat <- matrix(rnorm(300),100,3)
  tframe(dat) <- tframe(tspvector)
  ok <- is.tframed(dat)
  all.ok <- ok
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("tframe test 2 ... ")
  ok <- testEqual(tframe(dat), tframe(dat))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 3 ... ")
  ok <- all(c(1981,3) == start(tspvector))
  ok <- ok & all(c(1981,3) == start(dat))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 4 ... ")
  ok <- all(end(dat) == end(tspvector))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 5 ... ")
  ok <- Tobs(dat) == Tobs(tspvector)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 6 ... ")
  ok <- frequency(dat) == frequency(tspvector)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 7 ... ")
  #z <- tframed(dat, list(start=c(1961,2), frequency=12) )
  z <- ts(dat, start=c(1961,2), frequency=12) 
  ok <- is.tframed(z)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 8 ... ")
  z <- dat[10:90,]
  tframe(z) <- tfTruncate(tframe(dat), start=10, end=90)
  ok <- is.tframed(z)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 9 ... ")
  z <- tfTruncate(dat, start=10, end=90)
  ok <- is.tframed(z)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 10... ")
  #dat <- tframed(matrix(rnorm(300),100,3), list(start=c(1961,1), frequency=12))
  dat <- ts(matrix(rnorm(300),100,3), start=c(1961,1), frequency=12)
  z   <- tfwindow(dat, start=c(1963,2))
  zz  <- dat
  zz  <- tfwindow(zz, start=c(1963,2))
  zzz <- tfwindow(dat, start=c(1963,2))
  tframe(zzz) <- tframe(z)
  zzz  <- tframed(zzz, tframe(zzz))
  zzzz <- tframed(matrix(rnorm(300),100,3), tframe(dat))

  all.ok <- all.ok &  is.tframed(dat)
  all.ok <- all.ok &  is.tframed(z)
  all.ok <- all.ok &  is.tframed(zz)
  all.ok <- all.ok &  is.tframed(zzz)
  all.ok <- all.ok &  is.tframed(zzzz)
  all.ok <- all.ok &  all(z==zz)
  all.ok <- all.ok &  all(z==zzz)
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 11... ")
  ok <- all( time(dat) == time( tframed(dat, tframe(dat))))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 12... ")
#  z <- tsmatrix(1:10, 11:20)
  z <- tbind(1:10, 11:20)
  ok <-  all(start(z) ==1) & all( z== matrix(1:20, 10,2)) 
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 13... ")
  #dat <- tframed(matrix(rnorm(300),100,3), list(start=c(1961,1), frequency=12))
  dat <- ts(matrix(rnorm(300),100,3), start=c(1961,1), frequency=12)
  z <- tfwindow(dat, start=c(1963,2), end=c(1969,1))
  ok <-      all(start(dat)== earliestStart(dat, z))
  ok <- ok & all(    end(z) == earliestEnd  (dat, z))
  ok <- ok & all(start(z)   == latestStart  (dat, z))
  ok <- ok & all( end(dat) == latestEnd   (dat, z))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tframe test 14... ")
  dat <- tframed(matrix(rnorm(300),100,3), list(start=c(1961,1), frequency=12))
  z <- tfwindow(dat, start=c(1963,2), end=c(1969,1))
  ok <- testEqual(dat, splice(z, dat))
  ok <- ok & testEqual(tframe(dat), tframe(splice(z, dat)))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (synopsis) 
    {if (verbose) cat("All tframe tests completed")
     if (all.ok) cat(" OK\n") else cat(", some FAILED!\n") }
     
  if (all.ok) invisible(TRUE)  else stop("FAILED")
