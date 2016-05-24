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

  VARmodelB <- TSmodel(VARmodel)
  B <- t(chol(VARmodel$estimates$cov))
  VARmodelB$B <- array(B, c(1,dim(B)))  # has B != I
  VARmodelB <- setTSmodelParameters(VARmodelB)
  VARmodelB <- l(VARmodelB,VARmodel$data, warn=FALSE)

cat("dse test 8a ...\n")
  z  <- simulate(TSmodel(VARmodel), input=inputData(eg1.DSE.data.diff)) 
  zz <- simulate(TSmodel(VARmodel), rng=setRNG::getRNG(z), 
                     input=inputData(eg1.DSE.data.diff))
  ok <- testEqual(z, zz, fuzz=fuzz.small)

   if (!ok) 
     {
      all.ok <- FALSE  
     }

cat("dse test 8b ...\n")
  sigma <- solve(t(VARmodelB$model$B[1,,]) %*% VARmodelB$model$B[1,,])
  sigma <- (sigma + t(sigma))/2 # insure symetric - chol is sensitive
  zzz <- simulate(TSmodel(VARmodelB), rng=setRNG::getRNG(z), 
                     input=inputData(eg1.DSE.data.diff), Cov=sigma)
  ok <- testEqual(z, zzz, fuzz=fuzz.small)
  error <- max(abs(outputData(z) - outputData(zzz)))
   cat("max. error ", max(error), "\n")

   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error || !ok) 
     {printTestValue(outputData(z), digits=18)
      printTestValue(outputData(zzz), digits=18)
      all.ok <- FALSE  
     }

cat("dse test 8c ...\n")
  # next use estimates$cov
  z  <- simulate(VARmodel, input=inputData(eg1.DSE.data.diff)) 
  sigma <- VARmodel$estimates$cov
  sigma <- (sigma + t(sigma))/2 # insure symetric - chol is sensitive
  zz <- simulate(TSmodel(VARmodel), rng=setRNG::getRNG(z), 
                     input=inputData(eg1.DSE.data.diff), Cov=sigma)

  if ( ! testEqual(z, zz, fuzz=fuzz.small)) 
     {
      all.ok <- FALSE  
     }

  ok <- testEqual(summary(z)$estimates,
                   summary(zz)$estimates, fuzz=fuzz.small)

  if ( !ok) 
     {printTestValue(c(summary(z)$estimates),  digits=18)
      printTestValue(c(summary(zz)$estimates), digits=18)
      all.ok <- FALSE  
     }



  if (! all.ok) stop("some tests FAILED")

