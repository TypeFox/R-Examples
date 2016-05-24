require("dse")
Sys.info()
DSEversion()
 
fuzz <- 1e-6
digits <- 18
all.ok <- TRUE  

test.rng <- list(kind="Wichmann-Hill",seed=c(979,1479,1542),normal.kind="Box-Muller")

testModel <- SS(F=diag(1,3), H=matrix(c(1,0,0,1,0,0),2,3),
   Q=diag(0.5, 3, 3), R=diag(1.1, 2,2),
    description="test model", output.names=c("output 1", "output 2"))

z <- simulate(testModel, rng=test.rng)

if(! testEqual(z, simulate(testModel, rng=test.rng, compiled=FALSE),fuzz=1e-14))
     {cat("compiled and S versions of simulate differ!!!!!")
      all.ok <- FALSE  
     }

estModel <- estMaxLik(testModel, z)

#good <- 293.91258790365771 before simulate fix for w instead of e in non-innov models (Oct 2004)
good <-  340.556405433164741 
tst  <- estModel$estimates$like[1]
error <- max(abs(good - tst))
cat("max. error ", max(error), "\n")
 
if (any(is.na(error)) || any(is.nan(error)) || fuzz < error) 
     {printTestValue(c(tst), digits=18)
      all.ok <- FALSE  
     }
 
if (! all.ok) stop("some tests FAILED")

