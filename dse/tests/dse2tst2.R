  require("stats")
  require("dse") 
 #x11()
  dir <- tempdir()
  postscript(file=paste(dir,"/lite.out.ps", sep=""),  paper="letter",
      horizontal=FALSE, onefile=TRUE)
             # width=6, height=8, pointsize=10,
   Sys.info()
   DSEversion()
   setRNG::random.number.test() 




EvalEst.function.tests <- function(verbose=TRUE, synopsis=TRUE,
    fuzz.small=1e-14, fuzz.large=1e-8, graphics=TRUE)
{max.error <- NA
 data("eg1.DSE.data.diff", package="dse")
 
 if (synopsis & !verbose) cat("All EvalEst tests ...") 
 if (verbose) cat("EvalEst test 0 ... ")
  z <- eg1.DSE.data.diff
  z$input <- NULL
  mod1 <- TSmodel(estVARXar(z, re.add.means=FALSE, warn=FALSE))
  ok <- is.TSmodel(mod1)
  all.ok <- ok 
  if (verbose)  {if (ok) cat("ok\n") else  cat("failed!\n") }

  if (verbose) cat("EvalEst test 1 ... ")
  z <- estBlackBox1(eg1.DSE.data.diff, verbose=FALSE, max.lag=2)
  error <- max(abs(z$estimates$like[1]+4025.943051342767))
  ok <- is.TSestModel(z) &  (fuzz.large > error )
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  if (verbose) cat("EvalEst test 2 ... ")
  z <- estWtVariables(eg1.DSE.data.diff, c(1,10,10),
                        estimation="estVARXls")
  error <- max(abs(z$estimates$like[1]+4125.05572604540066)) 
  ok <- is.TSestModel(z) &   (fuzz.large > error)
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  if (verbose) cat("EvalEst test 3 ... ")
  z <- estSSMittnik(eg1.DSE.data.diff, max.lag=2, n=3)
  error <- max(abs(z$estimates$like[1]+3794.0394069904219))
  ok <- is.SS(z$model) &   (fuzz.large > error )
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  if (verbose) cat("EvalEst test 4 ... ")
  z <- l( MittnikReduction(z, criterion="taic", verbose=FALSE), 
         eg1.DSE.data.diff)
  error <- max(abs(z$estimates$like[1]+3795.6760513068380)) 
  ok <- is.SS(z$model)  &  (fuzz.large > error )
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  modSS<-z

  if (verbose) cat("EvalEst test 5 ... ")
  z <- featherForecasts( modSS,  from.periods=c(250,300))
  error <- max(abs
       (c(z$featherForecasts[[1]][286,],z$featherForecasts[[2]][336,])
       -c(-0.00092229286770808757701, -0.0086020067525247358164, 
           0.0043454851777852505565,  -0.0066741302949233430319,
          -0.0089398331205012854933,   0.0021769124280658046222)))
# 10* for change from svd to La.svd (but values not changed)
  ok <- 10*fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

  if (verbose) cat("EvalEst test 6 ... ")
  # previously end=c(1969,6) when .diff data had wrong start date
  outputData(modSS$data) <- tfwindow(outputData(modSS), end=c(1969,7))
  # it should be possible to do the following instead, but tsp seems to
  # sometimes get mixed up in forecast and cause System terminating: bad address
  # outputData(modSS$data) <- outputData(modSS$data)[1:100,]
  z <- forecast(modSS, percent=c(90,100,110))

# previously 136 below
  error <- max(abs(
    c(z$forecast[[1]][36,],z$forecast[[2]][36,], z$forecast[[3]][36,])
     -c(-0.00310702417651131587, -0.00604105559321206804,0.00214657444656118738,
      -0.00345224972784219028, -0.00671228396225603124,0.00238508249578931863,
      -0.00379747527917305948, -0.00738351233129999531,0.00262359054501745074)))
# 10* for change from svd to La.svd (but values not changed)
  ok <- 10*fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed! error = ", error,")\n") }

if (graphics) 
 {if (verbose) cat("EvalEst test 7 (graphics) ... ")
  ok <- EvalEst.graphics.tests(verbose=verbose, pause=TRUE)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }
 }

  if (synopsis) 
    {if (verbose) cat("All EvalEst tests completed")
     if (all.ok) cat(" OK\n")
     else    
       {cat(", some FAILED!")
        if(max.error > fuzz.small)
            cat(" max. error magnitude= ", max.error,")")
        cat("\n")
       }
    }

  if (all.ok) invisible(TRUE)  else stop("FAILED")
}

EvalEst.graphics.tests <- function(verbose=TRUE, synopsis=TRUE)
{# graphics tests do not do any value comparisons
  if (synopsis & !verbose) cat("EvalEst graphics tests ...")
  
 data("eg1.DSE.data.diff", package="dse")
 
 if (verbose) cat("  EvalEst graphics test 1 ...")

  # If no device is active then write to postscript file 
  if (dev.cur() == 1 )
      {postscript(file="zot.postscript.test.ps",width=6,height=6,pointsize=10,
                   onefile=FALSE, print.it=FALSE, append=FALSE)
       on.exit((function()
             {dev.off(); synchronize(1); rm("zot.postscript.test.ps")})())
      }

  data <- eg1.DSE.data.diff
  mod1 <- TSmodel(estVARXls(data,max.lag=3))
  modSS <- l(toSS(mod1),data)

  z <- featherForecasts( modSS,  from.periods=c(230,250))
  tfplot(z, start=c(1980,1))
  if (verbose) cat("ok\n")

  if (verbose) cat("  EvalEst graphics test 2 ...")
  z <- forecast(modSS, percent=c(90,100,110))
  tfplot(z, start=c(1985,1))
  if (verbose) cat("ok\n")

  if (synopsis) 
    {if (verbose) cat("All EvalEst graphics tests completed\n")
     else cat("completed\n")
    }
      

  invisible(TRUE)
}



   EvalEst.function.tests(verbose=TRUE, graphics=FALSE)  
   EvalEst.graphics.tests(verbose=TRUE)

unlink(dir, recursive=TRUE)
