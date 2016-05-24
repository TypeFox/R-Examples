 require("dse")
 Sys.info()
 DSEversion()
 data("eg1.DSE.data.diff", package="dse") 
 
 if (!is.TSdata(eg1.DSE.data.diff)) stop("Test data not found. Testing stopped.")
 
fuzz.small <- 1e-14
fuzz.large <- 1e-10
digits <- 18
all.ok <- TRUE  


test.rng <- list(kind="Wichmann-Hill",seed=c(979,1479,1542),normal.kind="Box-Muller")

VARmodel  <-  estVARXar(eg1.DSE.data.diff, re.add.means=FALSE, warn=FALSE)
SSmodel  <- toSS(VARmodel)


cat("dse test 2 ...\n")
   good <- VARmodel$estimates$like[1]
   tst  <- l(setArrays(SSmodel), eg1.DSE.data.diff,warn=FALSE)$estimates$like[1]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }

cat("dse test 3 ...\n")
   good <- VARmodel$estimates$like[1] 
   tst  <- l(setArrays(VARmodel), eg1.DSE.data.diff, warn=FALSE)$estimates$like[1]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }


cat("dse test 4 ...\n")
  ARMAmodel <- toARMA(SSmodel)
  good <- VARmodel$estimates$like[1]
   tst  <- l(ARMAmodel, eg1.DSE.data.diff, warn=FALSE)$estimates$like[1]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }


cat("dse test 5a...\n")
   good <- VARmodel$estimates$like[1]
   tst  <- l(ARMAmodel, eg1.DSE.data.diff,warn=FALSE)$estimates$like[1]
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.large < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }


cat("dse test 5b...\n")
  # longer input data is used by forecast in EvalEst
  data("egJofF.1dec93.data", package="dse")
  eg4.DSE.data<- egJofF.1dec93.data
  outputData(eg4.DSE.data) <- outputData(eg4.DSE.data, series=c(1,2,6,7))
  eg4.DSE.model <- estVARXls(eg4.DSE.data)
  longIn.data <- TSdata(
              input= ts(rbind(inputData(eg4.DSE.data), matrix(.1,10,1)), 
                       start=start(eg4.DSE.data),
                       frequency=frequency(eg4.DSE.data)),    
              output=outputData(eg4.DSE.data))
  seriesNames(longIn.data) <- seriesNames(eg4.DSE.data)
  z  <- l(TSmodel(eg4.DSE.model), longIn.data) 
  zz <- l(TSmodel(eg4.DSE.model), longIn.data, compiled=FALSE) 
  error <- max(abs(z$estimates$pred - zz$estimates$pred))
  #tfplot(z$estimates$pred, zz$estimates$pred)
  #z$estimates$pred[1:5,] ; zz$estimates$pred[1:5,]
  ok <-  fuzz.small > error 
  if (ok) cat("ok\n") else {
     max.error <- if (is.na(max.error)) error else max(error, max.error)
     cat("failed! error= ", error,"\n") 
     if(!testEqual(outputData(z), outputData(zz)))
         cat("output data comparison for l() and longIn.data failed.\n") 
     if(!testEqual(inputData(z), inputData(zz)))
         cat("input data comparison for l() and longIn.data failed.\n") 
     }
  all.ok <- all.ok & ok 


cat("dse test 6a...\n")
   good <- sort(Mod(roots(TSmodel(VARmodel),by.poly=TRUE)))
   tst  <- sort(Mod(roots(SSmodel)))
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }


cat("dse test 6b...\n")
   good <- sort(Mod(roots(SSmodel)))
   tst  <- sort(Mod(roots(TSmodel(VARmodel),by.poly=FALSE)))
   error <- max(abs(good-tst))
   cat("max. error ", max(error))

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }



  if (! all.ok) stop("some tests FAILED")

