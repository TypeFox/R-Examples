if(!require("stats"))   stop("this test requires stats.")
if(!require("dse"))  stop("this test requires dse.")
if(!require("EvalEst"))  stop("this test requires EvalEst.")
if(!require("setRNG"))stop("this test requires setRNG.")
 #x11()
  outfile <- tempfile("lite.out", tmpdir = tempdir(), fileext = ".ps")
  postscript(file=outfile,  paper="letter", horizontal=FALSE, onefile=TRUE)
             # width=6, height=8, pointsize=10,
   Sys.info()
   DSEversion()
   random.number.test() 


verbose <- TRUE 
synopsis <- TRUE
fuzz.small <- 1e-14 
fuzz.large <- 1e-8
graphics <- FALSE
max.error <- NA

 data("eg1.DSE.data.diff", package="dse")
# The seed is not important for most of these tests, but AIC eliminates all
#  parameters occassionally in some model selection tests.
 test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
# setRNG(test.rng)

  if (synopsis & !verbose) cat("All dse3 tests ...")
  if (verbose) cat("dse3 test 0 ... ")
  data <- eg1.DSE.data.diff
  inputData(data) <- NULL
  mod1 <- TSmodel(estVARXls(data))
  mod2 <- TSmodel(estVARXar(data, re.add.means=FALSE, warn=FALSE))
  ok <- is.TSmodel(mod1)
  all.ok <- ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 1 ... ")
  z <- MonteCarloSimulations(mod1, replications=5, quiet=TRUE)
  ok <- is.MonteCarloSimulations(z)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 2 ... ")
  ok <- testEqual(z, MonteCarloSimulations(mod1, replications=5,
                                     rng=getRNG(z), quiet=TRUE))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 3 ... ")
  z <- EstEval(mod1, replications=3,  estimation="estVARXls",
            estimation.args=NULL, criterion="TSmodel", quiet=TRUE)
  ok <- is.EstEval(z)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 4 ... ")
  zz <-summary(coef(z), verbose=FALSE)
  ok <- TRUE   # could be improved
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 5 ... ")
  zz <- summary(roots(z), verbose=FALSE)
  ok <- TRUE   # could be improved
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 6 ... ")
  zz <- coef(z)
  ok <- is.EstEval(zz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 7 ... ")
  zz <- roots(z)
  ok <- is.EstEval(zz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

#print(library.dynam())
  if (verbose) cat("dse3 test 8a... ")
  z <- horizonForecasts(mod1, data, horizons=c(6,12), discard.before=20)
  error <- max(abs( c(z$horizonForecasts[,100,])  -
 c(0.0048425425521641824594, 0.0031489473295282835973, 0.0037730234730729999594,
 0.0024354234760485438289, 0.0040593859721713481878, 0.0031982930612152113414)))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 8b... ")
  z <- horizonForecasts(l(toSS(mod1), data),
                         horizons=c(6,12), discard.before=20)
  error <- max(abs( c(z$horizonForecasts[,100,]) -
 c(0.0048425425521641824594, 0.0031489473295282844646, 0.0037730234730729995257,
 0.0024354234760485446963, 0.0040593859721713499225, 0.0031982930612152122088)))
  ok <- fuzz.small > error
  if (!ok) {if (is.na(max.error)) max.error <- error
            else max.error <- max(error, max.error)}
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("dse3 test 9 ... ")
  zzz<-l(mod1,simulate(mod1))
  zz<-forecastCov(zzz, discard.before=50, horizons=1:4)
  ok <- is.forecastCov(zz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 10... ")
  ok <- testEqual(zz, forecastCov(zzz$model, 
             data=zzz$data, discard.before=50, horizons=1:4))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 11... ")
  zz <-forecastCov(mod1,mod2, data=data, discard.before=30,
     zero=TRUE, trend=TRUE)

  ok <- is.forecastCov(zz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 12... ")
  zzz <-forecastCov(toSS(mod1),toSS(mod2), data=data, 
                 discard.before=30, zero=TRUE, trend=TRUE)

  ok <- testEqual(zz,zzz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 13... ")
  zz <-outOfSample.forecastCovEstimatorsWRTdata(data,
               estimation.methods = list(estVARXar= list(max.lag=2, warn=FALSE), 
                                         estVARXls= list(max.lag=2)))
  ok <- is.forecastCov(zz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 14... ")
  zz <- forecastCovWRTtrue(list(mod1,mod2),mod1, 
          pred.replications=2, quiet=TRUE, trend=NULL, zero=TRUE)
#          pred.replications=2, Spawn=FALSE, quiet=TRUE, trend=NULL, zero=TRUE)
  ok <- is.forecastCov(zz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 15... ")
  ok <- testEqual(zz, forecastCovWRTtrue(list(mod1,mod2),mod1, 
#          pred.replications=2, Spawn=if(exists(".SPAWN")) .SPAWN else FALSE,
          pred.replications=2,
	  quiet=TRUE, trend=NULL, zero=TRUE, rng=getRNG(zz)))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 16... ")
  zz <- forecastCovEstimatorsWRTtrue(mod1, 
#         Spawn=if(exists(".SPAWN")) .SPAWN else FALSE, quiet=TRUE, 
         quiet=TRUE, 
         estimation.methods=list(estVARXls=NULL, estVARXar=list(warn=FALSE)), 
         est.replications=2, pred.replications=2, rng=test.rng)
  ok <- is.forecastCov(zz)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse3 test 17... ")
# Next seems to cause a problem in Splus when .SPAWN is TRUE above although it may
#  work with default rng (at least it used to) but test.rng is now set
#  to give same results as in R.
#  ok <- testEqual(zz, forecastCovEstimatorsWRTtrue(mod1, Spawn=FALSE, 
  ok <- testEqual(zz, forecastCovEstimatorsWRTtrue(mod1,  
           estimation.methods=list(estVARXls=NULL,estVARXar=list(warn=FALSE)), 
           est.replications=2, pred.replications=2, rng=getRNG(zz)))
  all.ok <- all.ok & ok 
  if (verbose) if (ok) cat("ok\n") else cat("failed!\n") 

  if (graphics)
      {ok <- dse3.graphics.tests(verbose=verbose,  pause=FALSE)
       all.ok <- all.ok & ok 
       if (verbose) cat("dse3 test 18 (graphics) ... ")
       if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }
      }

  if (synopsis) 
    {if (verbose) cat("All dse3 tests completed")
     if (all.ok) cat(" OK\n") 
     else  cat(" some FAILED! max.error = ", max.error,"\n")
    }

  if (all.ok) invisible(TRUE)  else stop("FAILED")



# dse3.graphics.tests 

verbose <- TRUE
synopsis <- TRUE

data("eg1.DSE.data.diff", package="dse")
 
  if (synopsis & !verbose) cat("dse3 graphics tests ...")
  if (verbose) cat("  dse3 graphics test 1 ...")
  # If no device is active then write to postscript file 
  if ( dev.cur() == 1 )
      {postscript(file="zot.postscript.test.ps",
                   width=6,height=6,pointsize=10,
                   onefile=FALSE, print.it=FALSE, append=FALSE)
       on.exit((function()
            {dev.off(); synchronize(1); rm("zot.postscript.test.ps")})())
      }

# The seed is not important for most of these tests, but AIC eliminates all
#  parameters occassionally in some model selection tests.
 test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))

  data <- eg1.DSE.data.diff
  inputData(data) <- NULL
  outputData(data) <- outputData(data, series=1)  # [,1,drop=FALSE]
  mod1 <- TSmodel(estVARXls(data,max.lag=3))
  mod2 <- TSmodel(estVARXar(data,max.lag=3, aic=FALSE, warn=FALSE))

  z <- EstEval(mod1, replications=10,  estimation="estVARXls",
            estimation.args=list(max.lag=3), criterion="TSmodel", quiet=TRUE)
  distribution(coef(z)) 
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse3 graphics test 2 ...")
  distribution(roots(z))
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse3 graphics test 3 ...")
  z <- horizonForecasts(mod1, data, horizons=c(6,12), discard.before=20)
  tfplot(z, start=c(1985,1))
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse3 graphics test 4 ...")
  zz <-forecastCov(mod1,mod2, data=data,
                     discard.before=10, zero=TRUE, trend=TRUE)
  tfplot(zz)
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse3 graphics test 5 ...")
  tfplot(zz, select.cov=c(1), select.trend=FALSE)
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse3 graphics test 6 ...")

  data <- eg1.DSE.data.diff
  inputData(data) <- NULL
# next causes Error ... all lags eliminated by AIC order selection.
#  outputData(data) <- outputData(data, series=1)  # [,1,drop=FALSE]
  zz <-outOfSample.forecastCovEstimatorsWRTdata(data,
               estimation.methods = list(estVARXar=list(max.lag=2,warn=FALSE),
                                         estVARXls=list(max.lag=2))) 
  tfplot(zz, series=c(1))
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse3 graphics test 7 ...")
  zz <- forecastCovWRTtrue(list(mod1,mod2),mod1, rng=test.rng,
               pred.replications=2, 
#	       Spawn=if(exists(".SPAWN")) .SPAWN else FALSE, 
	       trend=NULL, zero=TRUE, quiet=TRUE)
  tfplot(zz, select.cov=c(1))
  if (verbose) cat("ok\n")

  if (verbose) cat("  dse3 graphics test 8 ...")
  zz <- forecastCovEstimatorsWRTtrue(mod1, 
#         Spawn=if(exists(".SPAWN")) .SPAWN else FALSE,  
	 rng=test.rng,
         estimation.methods=list(estVARXls=NULL,estVARXls=list(max.lag=2)), 
#        estimation.methods=list(estVARXls=NULL,estVARXar=NULL), 
         est.replications=2, pred.replications=2, quiet=TRUE)
  tfplot(zz, select.cov=c(1))
  if (verbose) cat("ok\n")

  if (synopsis) 
    {if (verbose) cat("All dse3 graphics tests completed\n")
     else cat("completed\n")
    }
      
unlink(outfile)
