if(!require("stats"))   stop("this test requires stats.")
if(!require("EvalEst"))  stop("this test requires EvalEst.")
if(!require("setRNG"))stop("this test requires setRNG.")
 #x11()
  outfile <- tempfile("lite.out", tmpdir = tempdir(), fileext = ".ps")
  postscript(file=outfile,  paper="letter", horizontal=FALSE, onefile=TRUE)
             # width=6, height=8, pointsize=10,
   Sys.info()
   DSEversion()
   random.number.test() 




dse4.function.tests <- function(verbose=TRUE, synopsis=TRUE, 
		fuzz.small=1e-14, fuzz.large=1e-7, graphics=TRUE)
{max.error <- 0
 data("eg1.DSE.data.diff", package="dse")
 
 if (synopsis & !verbose) cat("All dse4 tests ...") 
 if (verbose) cat("dse4 test 1 ... ")
  z <- stripMine(eg1.DSE.data.diff, essential.data=c(1,2),
                   estimation.methods=list(estVARXls=list(max.lag=3)))
  ok <- is.forecastCovEstimatorsWRTdata.subsets(z)
  all.ok <-  ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse4 test 2 ... ")
  z1 <- z$multi.model[[
       selectForecastCov(z, select.cov.best=1, verbose=FALSE)$selection.index[2]]]
  subdata <- TSdata(output=outputData(eg1.DSE.data.diff, series=1:3))
  z2 <- estimateModels(subdata, estimation.sample =182, quiet = TRUE, 
           estimation.methods = list(estVARXls=list(max.lag=3)))
  outputData(subdata) <- outputData(subdata)[1:182,,drop=FALSE]
#  inputData(subdata)  <- inputData(subdata) [1:182,,drop=FALSE] not in subdata
  z3 <- estVARXls(subdata, max.lag=3)
  # Had to add fuzz=1e-15 in next two, default =0 no longer worked April 2012 
  # with R-2.15.0 on 32 bit Ubunti (with new BLAS recently installed).
  ok <-      testEqual(z2$multi.model[[1]],z3$model, fuzz=1e-15)
  ok <- ok & testEqual(z2$multi.model[[1]],  z1, fuzz=1e-15)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("dse4 test 3 ... ")
  cat("skipping test 4 (requires stepwise).")
  if (FALSE) 
 {
   all.data <- TSdata(input=eg1.DSE.data.diff$output, 
                   output=eg1.DSE.data.diff$input )
   umodel <- build.input.models(all.data, max.lag=2)
   umodel <- build.diagonal.model(umodel)
   z  <- TSdata(output=outputData(all.data), 
                input=inputData(all.data, series=1:2))
	# previously ??input=inputData(extract(all.data, outputs=1, inputs=1:2)))
  ymodel <- estVARXls(z, max.lag=3)$model 
  z <- ymodel$C
  ymodel$C <- array(0, c(dim(z)[1:2], nseriesOutput(umodel))) 
  ymodel$C[1:(dim(z)[1]), 1:(dim(z)[2]), 1:(dim(z)[3])] <- z 
  sim.data <- genMineData(umodel, ymodel,
    rng= list(kind="default",seed=c(21,46,16,12, 51, 2, 31, 8, 42, 60, 7, 3)) )
  m.step <- mineStepwise(sim.data, method="backward", plot.=FALSE)
  error <- max(abs(m.step$stepwise$rss[c(1,27)] -
                c(5.65168201726030333e+01, 3.06441399376576440e+06)))
  # previously  c(47.537312899054931847,   4088283.2706551752053)))
  ok <- fuzz.large > error
  if (!ok) max.error <- max(error, max.error)
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }
 }

  if (graphics)
   {ok <- dse4.graphics.tests(verbose=verbose, pause=FALSE)
    all.ok <- all.ok & ok 
   }

  if (synopsis) 
    {if (verbose) cat("All dse4 tests completed")
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

dse4.graphics.tests <- function(verbose=TRUE, synopsis=TRUE)
{ if (synopsis & !verbose) cat("dse4 graphics tests ...")
  if (verbose) cat("  dse4 graphics test 1 ...")

  # If no device is active then write to postscript file 
  if ( dev.cur() == 1 )
      {postscript(file="zot.postscript.test.ps",width=6,height=6,pointsize=10,
                   onefile=FALSE, print.it=FALSE, append=FALSE)
       on.exit((function()
             {dev.off(); synchronize(1); rm("zot.postscript.test.ps")})())
      }

  z <- stripMine(eg1.DSE.data.diff, essential.data=c(1,2),
                   estimation.methods=list(estVARXls=list(max.lag=3)))
  zz <- tfplot(z)

  if (verbose) cat("ok\n")

  if (synopsis) 
    {if (verbose) cat("All dse4 graphics tests completed\n")
     else cat("completed\n")
    }
      
  invisible(TRUE)
}



   dse4.function.tests(verbose=TRUE, graphics=FALSE) 
   dse4.graphics.tests(verbose=TRUE)  #     test 3 needs stepwise

unlink(outfile)
