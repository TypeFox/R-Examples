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

cat("dse test 9 ...\n")
  z  <- simulate(SSmodel, input=inputData(eg1.DSE.data.diff)) 
  ok <- testEqual(z,simulate(SSmodel, rng=setRNG::getRNG(z), 
                      input=inputData(eg1.DSE.data.diff)))
  if (!ok) {all.ok <- FALSE ; cat(ok, "\n")}

  ok <- testEqual(summary(z)$estimates,
                   summary(z)$estimates, fuzz=fuzz.small)
  if (!ok) {all.ok <- FALSE ; cat(ok, "\n")}



cat("dse test 10...\n")

  ok <- stability(SSmodel, verbose=FALSE)
  if (!ok) {all.ok <- FALSE ; cat(ok, "\n")}


cat("dse test 11...\n")

   scale.fac <- diag(1:3)
   scale.fac[1,3] <-.5
   scale.pred <- VARmodel$estimates$pred %*% t(scale.fac)
   scale.fac <- list(output=scale.fac)

   good <- scale.pred
   tst  <- l(scale(VARmodel$model, scale=scale.fac), 
          scale(eg1.DSE.data.diff, scale=scale.fac), warn=FALSE)$estimates$pred
   error <- max(abs(good - tst))
   cat("max. error ", max(error), "\n")
 
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }


cat("dse test 12...\n")

   good <- scale.pred
   tst  <- l(scale(SSmodel, scale=scale.fac), 
             scale(eg1.DSE.data.diff, scale=scale.fac))$estimates$pred
   error <- max(abs(good - tst))
   cat("max. error ", max(error), "\n")
 
   if (any(is.na(error)) || any(is.nan(error)) || fuzz.small < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }


cat("dse test 13...\n")

  z <- eg1.DSE.data.diff
  ok <- testEqual(z,
      TSdata(output=outputData(combine(z,z), series=seq(nseriesOutput(z))),
              input= inputData(combine(z,z), series=seq( nseriesInput(z))))) 
 
  if (!ok) {all.ok <- FALSE ; cat(ok, "\n")}


  if (! all.ok) stop("some tests FAILED")

