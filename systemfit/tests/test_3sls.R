library( systemfit )
options( digits = 3 )

data( "Kmenta" )
useMatrix <- FALSE

demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
inst <- ~ income + farmPrice + trend
inst1  <- ~ income + farmPrice
instlist <- list( inst1, inst )
system <- list( demand = demand, supply = supply )
restrm <- matrix(0,1,7)  # restriction matrix "R"
restrm[1,3] <-  1
restrm[1,7] <- -1
restrict <- "demand_income - supply_trend = 0"
restr2m <- matrix(0,2,7)  # restriction matrix "R" 2
restr2m[1,3] <-  1
restr2m[1,7] <- -1
restr2m[2,2] <- -1
restr2m[2,5] <-  1
restr2q <- c( 0, 0.5 )  # restriction vector "q" 2
restrict2 <- c( "demand_income - supply_trend = 0",
   "- demand_price + supply_price = 0.5" )
tc <- matrix(0,7,6)
tc[1,1] <- 1
tc[2,2] <- 1
tc[3,3] <- 1
tc[4,4] <- 1
tc[5,5] <- 1
tc[6,6] <- 1
tc[7,3] <- 1
restr3m <- matrix(0,1,6)  # restriction matrix "R" 2
restr3m[1,2] <- -1
restr3m[1,5] <-  1
restr3q <- c( 0.5 )  # restriction vector "q" 2
restrict3 <- "- C2 + C5 = 0.5"


## *************** 3SLS estimation ************************
fit3sls <- list()
formulas <- c( "GLS", "IV", "Schmidt", "GMM", "EViews" )
for( i in seq( along = formulas ) ) {
   fit3sls[[ i ]] <- list()

   print( "***************************************************" )
   print( paste( "3SLS formula:", formulas[ i ] ) )
   print( "************* 3SLS *********************************" )
   fit3sls[[ i ]]$e1 <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, method3sls = formulas[ i ], useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e1 ) )

   print( "********************* 3SLS EViews-like *****************" )
   fit3sls[[ i ]]$e1e <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "noDfCor", method3sls = formulas[ i ],
      useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e1e, useDfSys = TRUE ) )

   print( "********************* 3SLS with methodResidCov = Theil *****************" )
   fit3sls[[ i ]]$e1c <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "Theil", method3sls = formulas[ i ],
      x = TRUE, useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e1c, useDfSys = TRUE ) )

   print( "*************** W3SLS with methodResidCov = Theil *****************" )
   fit3sls[[ i ]]$e1wc <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "Theil", method3sls = formulas[ i ],
      residCovWeighted = TRUE, useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e1wc, useDfSys = TRUE ) )


   print( "*************** 3SLS with restriction *****************" )
   fit3sls[[ i ]]$e2 <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.matrix = restrm, method3sls = formulas[ i ],
      x = TRUE, useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e2 ) )
   # the same with symbolically specified restrictions
   fit3sls[[ i ]]$e2Sym <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.matrix = restrict, method3sls = formulas[ i ],
      x = TRUE, useMatrix = useMatrix )
   print( all.equal( fit3sls[[ i ]]$e2, fit3sls[[ i ]]$e2Sym ) )

   print( "************** 3SLS with restriction (EViews-like) *****************" )
   fit3sls[[ i ]]$e2e <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "noDfCor", restrict.matrix = restrm,
      method3sls = formulas[ i ], useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e2e, useDfSys = TRUE ) )
   print( nobs( fit3sls[[i]]$e2e ))

   print( "*************** W3SLS with restriction *****************" )
   fit3sls[[ i ]]$e2w <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.matrix = restrm, method3sls = formulas[ i ],
      residCovWeighted = TRUE, useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e2w ) )


   print( "*************** 3SLS with restriction via restrict.regMat ********************" )
   fit3sls[[ i ]]$e3 <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.regMat = tc, method3sls = formulas[ i ],
      useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e3 ) )

   print( "*************** 3SLS with restriction via restrict.regMat (EViews-like) *******" )
   fit3sls[[ i ]]$e3e <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "noDfCor", restrict.regMat = tc,
      method3sls = formulas[ i ], x = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e3e, useDfSys = TRUE ) )

   print( "**** W3SLS with restriction via restrict.regMat (EViews-like) ****" )
   fit3sls[[ i ]]$e3we <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "noDfCor", restrict.regMat = tc,
      method3sls = formulas[ i ], residCovWeighted = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e3we, useDfSys = TRUE ) )


   print( "*************** 3SLS with 2 restrictions **********************" )
   fit3sls[[ i ]]$e4 <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.matrix = restr2m, restrict.rhs = restr2q,
      method3sls = formulas[ i ], x = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e4 ) )
   # the same with symbolically specified restrictions
   fit3sls[[ i ]]$e4Sym <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.matrix = restrict2, method3sls = formulas[ i ],
      x = TRUE, useMatrix = useMatrix )
   print( all.equal( fit3sls[[ i ]]$e4, fit3sls[[ i ]]$e4Sym ) )

   print( "*************** 3SLS with 2 restrictions (EViews-like) ************" )
   fit3sls[[ i ]]$e4e <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "noDfCor", restrict.matrix = restr2m,
      restrict.rhs = restr2q, method3sls = formulas[ i ],
      useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e4e, useDfSys = TRUE ) )

   print( "********** W3SLS with 2 (symbolic) restrictions ***************" )
   fit3sls[[ i ]]$e4wSym <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.matrix = restrict2, method3sls = formulas[ i ],
      residCovWeighted = TRUE, useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e4wSym ) )


   print( "*************** 3SLS with 2 restrictions via R and restrict.regMat **********" )
   fit3sls[[ i ]]$e5 <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.regMat = tc, restrict.matrix = restr3m,
      restrict.rhs = restr3q, method3sls = formulas[ i ],
      useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e5 ) )
   # the same with symbolically specified restrictions
   fit3sls[[ i ]]$e5Sym <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.regMat = tc, restrict.matrix = restrict3,
      method3sls = formulas[ i ], useMatrix = useMatrix )
   print( all.equal( fit3sls[[ i ]]$e5, fit3sls[[ i ]]$e5Sym ) )

   print( "******** 3SLS with 2 restrictions via R and restrict.regMat (EViews-like)*****" )
   fit3sls[[ i ]]$e5e <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.regMat = tc, methodResidCov = "noDfCor",
      restrict.matrix = restr3m, restrict.rhs = restr3q,
      method3sls = formulas[ i ], x = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e5e, useDfSys = TRUE ) )

   print( "*** W3SLS with 2 restrictions via R and restrict.regMat (EViews-like) ***" )
   fit3sls[[ i ]]$e5we <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.regMat = tc, methodResidCov = "noDfCor",
      restrict.matrix = restr3m, restrict.rhs = restr3q, method3sls = formulas[ i ],
      residCovWeighted = TRUE, useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$e5we, useDfSys = TRUE ) )

   ## *********** estimations with a single regressor ************
   fit3sls[[ i ]]$S1 <- systemfit(
      list( farmPrice ~ consump - 1, price ~ consump + trend ), "3SLS",
      data = Kmenta, inst = ~ trend + income, useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$S1 ) )
   fit3sls[[ i ]]$S2 <- systemfit(
      list( consump ~ farmPrice - 1, consump ~ trend - 1 ), "3SLS",
      data = Kmenta, inst = ~ price + income, useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$S2 ) )
   fit3sls[[ i ]]$S3 <- systemfit(
      list( consump ~ trend - 1, farmPrice ~ trend - 1 ), "3SLS",
      data = Kmenta, inst = instlist, useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$S3 ) )
   fit3sls[[ i ]]$S4 <- systemfit(
      list( consump ~ farmPrice - 1, price ~ trend - 1 ), "3SLS",
      data = Kmenta, inst = ~ farmPrice + trend + income,
      restrict.matrix = matrix( c( 1, -1 ), nrow = 1 ), useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$S4 ) )
   fit3sls[[ i ]]$S5 <- systemfit(
      list( consump ~ 1, price ~ 1 ), "3SLS",
      data = Kmenta, inst = ~ income, useMatrix = useMatrix )
   print( summary( fit3sls[[ i ]]$S5 ) )
}

## ******************** iterated 3SLS **********************
fit3slsi <- list()
formulas <- c( "GLS", "IV", "Schmidt", "GMM", "EViews" )
for( i in seq( along = formulas ) ) {
   fit3slsi[[ i ]] <- list()

   print( "***************************************************" )
   print( paste( "3SLS formula:", formulas[ i ] ) )
   print( "************* 3SLS *********************************" )
   fit3slsi[[ i ]]$e1 <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, method3sls = formulas[ i ], maxiter = 100,
      useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e1 ) )

   print( "********************* iterated 3SLS EViews-like ****************" )
   fit3slsi[[ i ]]$e1e <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "noDfCor", method3sls = formulas[ i ],
      maxiter = 100, useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e1e, useDfSys = TRUE ) )

   print( "************** iterated 3SLS with methodResidCov = Theil **************" )
   fit3slsi[[ i ]]$e1c <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "Theil", method3sls = formulas[ i ],
      maxiter = 100, x = TRUE, useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e1c, useDfSys = TRUE ) )

   print( "**************** iterated W3SLS EViews-like ****************" )
   fit3slsi[[ i ]]$e1we <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "noDfCor", method3sls = formulas[ i ],
      maxiter = 100, residCovWeighted = TRUE, useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e1we, useDfSys = TRUE ) )


   print( "******* iterated 3SLS with restriction *****************" )
   fit3slsi[[ i ]]$e2 <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.matrix = restrm, method3sls = formulas[ i ],
      maxiter = 100, x = TRUE, useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e2 ) )

   print( "********* iterated 3SLS with restriction (EViews-like) *********" )
   fit3slsi[[ i ]]$e2e <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "noDfCor", restrict.matrix = restrm,
      method3sls = formulas[ i ], maxiter = 100, useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e2e, useDfSys = TRUE ) )

   print( "******** iterated W3SLS with restriction (EViews-like) *********" )
   fit3slsi[[ i ]]$e2we <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "noDfCor", restrict.matrix = restrm,
      method3sls = formulas[ i ], maxiter = 100, residCovWeighted = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e2we, useDfSys = TRUE ) )


   print( "********* iterated 3SLS with restriction via restrict.regMat *****************" )
   fit3slsi[[ i ]]$e3 <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.regMat = tc, method3sls = formulas[ i ],
      maxiter = 100, useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e3 ) )

   print( "********* iterated 3SLS with restriction via restrict.regMat (EViews-like) ***" )
   fit3slsi[[ i ]]$e3e <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "noDfCor", restrict.regMat = tc,
      method3sls = formulas[ i ], maxiter = 100, x = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e3e, useDfSys = TRUE ) )

   print( "***** iterated W3SLS with restriction via restrict.regMat ********" )
   fit3slsi[[ i ]]$e3w <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.regMat = tc, method3sls = formulas[ i ], maxiter = 100,
      residCovWeighted = TRUE, x = TRUE, useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e3w ) )


   print( "******** iterated 3SLS with 2 restrictions *********************" )
   fit3slsi[[ i ]]$e4 <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.matrix = restr2m, restrict.rhs = restr2q,
      method3sls = formulas[ i ], maxiter = 100, x = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e4 ) )

   print( "********* iterated 3SLS with 2 restrictions (EViews-like) *******" )
   fit3slsi[[ i ]]$e4e <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "noDfCor", restrict.matrix = restr2m,
      restrict.rhs = restr2q, method3sls = formulas[ i ], maxiter = 100,
      useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e4e, useDfSys = TRUE ) )

   print( "******** iterated W3SLS with 2 restrictions (EViews-like) *******" )
   fit3slsi[[ i ]]$e4we <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, methodResidCov = "noDfCor", restrict.matrix = restr2m,
      restrict.rhs = restr2q, method3sls = formulas[ i ], maxiter = 100,
      residCovWeighted = TRUE, useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e4we, useDfSys = TRUE ) )


   print( "******** iterated 3SLS with 2 restrictions via R and restrict.regMat *********" )
   fit3slsi[[ i ]]$e5 <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.regMat = tc, restrict.matrix = restr3m,
      restrict.rhs = restr3q, method3sls = formulas[ i ], maxiter = 100,
      useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e5 ) )

   print( "*** iterated 3SLS with 2 restrictions via R and restrict.regMat (EViews-like)**" )
   fit3slsi[[ i ]]$e5e <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.regMat = tc, methodResidCov = "noDfCor",
      restrict.matrix = restr3m, restrict.rhs = restr3q,
      method3sls = formulas[ i ], maxiter = 100, useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e5e, useDfSys = TRUE ) )

   print( "** iterated W3SLS with 2 restrictions via R and restrict.regMat ***" )
   fit3slsi[[ i ]]$e5w <- systemfit( system, "3SLS", data = Kmenta,
      inst = inst, restrict.regMat = tc, restrict.matrix = restr3m,
      restrict.rhs = restr3q, method3sls = formulas[ i ], maxiter = 100,
      residCovWeighted = TRUE, x = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3slsi[[ i ]]$e5w ) )
}

## **************** 3SLS with different instruments *************
fit3slsd <- list()
formulas <- c( "GLS", "IV", "Schmidt", "GMM", "EViews" )
for( i in seq( along = formulas ) ) {
   fit3slsd[[ i ]] <- list()

   print( "***************************************************" )
   print( paste( "3SLS formula:", formulas[ i ] ) )
   print( "************* 3SLS with different instruments **************" )
   fit3slsd[[ i ]]$e1 <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, method3sls = formulas[ i ], useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e1 ) )

   print( "******* 3SLS with different instruments (EViews-like) **********" )
   fit3slsd[[ i ]]$e1e <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, methodResidCov = "noDfCor", method3sls = formulas[ i ],
      useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e1e, useDfSys = TRUE ) )

   print( "**** 3SLS with different instruments and methodResidCov = Theil ***" )
   fit3slsd[[ i ]]$e1c <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, methodResidCov = "Theil", method3sls = formulas[ i ],
      x = TRUE, useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e1c, useDfSys = TRUE ) )

   print( "************* W3SLS with different instruments **************" )
   fit3slsd[[ i ]]$e1w <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, method3sls = formulas[ i ], residCovWeighted = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e1w ) )


   print( "******* 3SLS with different instruments and restriction ********" )
   fit3slsd[[ i ]]$e2 <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, restrict.matrix = restrm, method3sls = formulas[ i ],
      x = TRUE, useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e2 ) )

   print( "** 3SLS with different instruments and restriction (EViews-like) *" )
   fit3slsd[[ i ]]$e2e <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, methodResidCov = "noDfCor", restrict.matrix = restrm,
      method3sls = formulas[ i ], useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e2e, useDfSys = TRUE ) )

   print( "** W3SLS with different instruments and restriction (EViews-like) *" )
   fit3slsd[[ i ]]$e2we <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, methodResidCov = "noDfCor", restrict.matrix = restrm,
      method3sls = formulas[ i ], residCovWeighted = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e2we, useDfSys = TRUE ) )


   print( "** 3SLS with different instruments and restriction via restrict.regMat *******" )
   fit3slsd[[ i ]]$e3 <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, restrict.regMat = tc, method3sls = formulas[ i ],
      useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e3 ) )

   print( "3SLS with different instruments with restriction via restrict.regMat (EViews-like)" )
   fit3slsd[[ i ]]$e3e <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, methodResidCov = "noDfCor", restrict.regMat = tc,
      method3sls = formulas[ i ], x = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e3e, useDfSys = TRUE ) )

   print( "** W3SLS with different instr. and restr. via restrict.regMat ****" )
   fit3slsd[[ i ]]$e3w <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, restrict.regMat = tc, method3sls = formulas[ i ],
      residCovWeighted = TRUE, x = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e3w ) )


   print( "****** 3SLS with different instruments and 2 restrictions *********" )
   fit3slsd[[ i ]]$e4 <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, restrict.matrix = restr2m, restrict.rhs = restr2q,
      method3sls = formulas[ i ], x = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e4 ) )

   print( "** 3SLS with different instruments and 2 restrictions (EViews-like) *" )
   fit3slsd[[ i ]]$e4e <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, methodResidCov = "noDfCor", restrict.matrix = restr2m,
      restrict.rhs = restr2q, method3sls = formulas[ i ],
      useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e4e, useDfSys = TRUE ) )

   print( "**** W3SLS with different instruments and 2 restrictions *********" )
   fit3slsd[[ i ]]$e4w <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, restrict.matrix = restr2m, restrict.rhs = restr2q,
      method3sls = formulas[ i ], residCovWeighted = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e4w ) )


   print( " 3SLS with different instruments with 2 restrictions via R and restrict.regMat" )
   fit3slsd[[ i ]]$e5 <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, restrict.regMat = tc, restrict.matrix = restr3m,
      restrict.rhs = restr3q, method3sls = formulas[ i ],
      useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e5 ) )

   print( "3SLS with diff. instruments and 2 restr. via R and restrict.regMat (EViews-like)" )
   fit3slsd[[ i ]]$e5e <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, restrict.regMat = tc, methodResidCov = "noDfCor",
      restrict.matrix = restr3m, restrict.rhs = restr3q,
      method3sls = formulas[ i ], x = TRUE,
      useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e5e, useDfSys = TRUE ) )

   print( "W3SLS with diff. instr. and 2 restr. via R and restrict.regMat (EViews-like)" )
   fit3slsd[[ i ]]$e5we <- systemfit( system, "3SLS", data = Kmenta,
      inst = instlist, restrict.regMat = tc, methodResidCov = "noDfCor",
      restrict.matrix = restr3m, restrict.rhs = restr3q, method3sls = formulas[ i ],
      residCovWeighted = TRUE, useMatrix = useMatrix )
   print( summary( fit3slsd[[ i ]]$e5we, useDfSys = TRUE ) )
}


## **************** shorter summaries **********************
print( summary( fit3sls[[ 2 ]]$e1c, equations = FALSE ) )

print( summary( fit3sls[[ 3 ]]$e2e ), residCov = FALSE )

print( summary( fit3sls[[ 4 ]]$e3, useDfSys = FALSE ), residCov = FALSE )

print( summary( fit3sls[[ 5 ]]$e4e, equations = FALSE ),
   equations = FALSE )

print( summary( fit3sls[[ 1 ]]$e4wSym, residCov = FALSE ),
   equations = FALSE )

print( summary( fit3sls[[ 2 ]]$e5, residCov = FALSE ), residCov = TRUE )

print( summary( fit3slsi[[ 3 ]]$e3e, residCov = FALSE,
   equations = FALSE ) )

print( summary( fit3slsi[[ 4 ]]$e1we ), equations = FALSE, residCov = FALSE )

print( summary( fit3slsd[[ 5 ]]$e4, residCov = FALSE ) )

print( summary( fit3slsd[[ 1 ]]$e2we, equations = FALSE ) )


## ****************** residuals **************************
print( residuals( fit3sls[[ 1 ]]$e1c ) )
print( residuals( fit3sls[[ 1 ]]$e1c$eq[[ 1 ]] ) )

print( residuals( fit3sls[[ 4 ]]$e1wc ) )
print( residuals( fit3sls[[ 4 ]]$e1wc$eq[[ 1 ]] ) )

print( residuals( fit3sls[[ 2 ]]$e2e ) )
print( residuals( fit3sls[[ 2 ]]$e2e$eq[[ 2 ]] ) )

print( residuals( fit3sls[[ 3 ]]$e3 ) )
print( residuals( fit3sls[[ 3 ]]$e3$eq[[ 1 ]] ) )

print( residuals( fit3sls[[ 4 ]]$e4e ) )
print( residuals( fit3sls[[ 4 ]]$e4e$eq[[ 2 ]] ) )

print( residuals( fit3sls[[ 5 ]]$e5 ) )
print( residuals( fit3sls[[ 5 ]]$e5$eq[[ 1 ]] ) )

print( residuals( fit3slsi[[ 2 ]]$e3e ) )
print( residuals( fit3slsi[[ 2 ]]$e3e$eq[[ 1 ]] ) )

print( residuals( fit3slsi[[ 1 ]]$e2we ) )
print( residuals( fit3slsi[[ 1 ]]$e2we$eq[[ 1 ]] ) )

print( residuals( fit3slsd[[ 3 ]]$e4 ) )
print( residuals( fit3slsd[[ 3 ]]$e4$eq[[ 2 ]] ) )

print( residuals( fit3slsd[[ 5 ]]$e5we ) )
print( residuals( fit3slsd[[ 5 ]]$e5we$eq[[ 2 ]] ) )


## *************** coefficients *********************
print( round( coef( fit3sls[[ 3 ]]$e1c ), digits = 6 ) )
print( round( coef( fit3sls[[ 4 ]]$e1c$eq[[ 2 ]] ), digits = 6 ) )

print( round( coef( fit3slsi[[ 4 ]]$e2 ), digits = 6 ) )
print( round( coef( fit3slsi[[ 5 ]]$e2$eq[[ 1 ]] ), digits = 6 ) )

print( round( coef( fit3sls[[ 2 ]]$e2w ), digits = 6 ) )
print( round( coef( fit3sls[[ 3 ]]$e2w$eq[[ 1 ]] ), digits = 6 ) )

print( round( coef( fit3slsd[[ 5 ]]$e3e ), digits = 6 ) )
print( round( coef( fit3slsd[[ 5 ]]$e3e, modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( fit3slsd[[ 1 ]]$e3e$eq[[ 2 ]] ), digits = 6 ) )

print( round( coef( fit3slsd[[ 4 ]]$e3w ), digits = 6 ) )
print( round( coef( fit3slsd[[ 4 ]]$e3w, modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( fit3slsd[[ 5 ]]$e3w$eq[[ 2 ]] ), digits = 6 ) )

print( round( coef( fit3sls[[ 1 ]]$e4 ), digits = 6 ) )
print( round( coef( fit3sls[[ 2 ]]$e4$eq[[ 1 ]] ), digits = 6 ) )

print( round( coef( fit3slsi[[ 2 ]]$e4we ), digits = 6 ) )
print( round( coef( fit3slsi[[ 1 ]]$e4we$eq[[ 1 ]] ), digits = 6 ) )

print( round( coef( fit3slsi[[ 2 ]]$e5e ), digits = 6 ) )
print( round( coef( fit3slsi[[ 2 ]]$e5e, modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( fit3slsi[[ 3 ]]$e5e$eq[[ 2 ]] ), digits = 6 ) )


## *************** coefficients with stats *********************
print( round( coef( summary( fit3sls[[ 3 ]]$e1c, useDfSys = FALSE ) ),
   digits = 6 ) )
print( round( coef( summary( fit3sls[[ 4 ]]$e1c$eq[[ 2 ]], useDfSys = FALSE ) ),
   digits = 6 ) )

print( round( coef( summary( fit3slsd[[ 2 ]]$e1w, useDfSys = FALSE ) ),
   digits = 6 ) )
print( round( coef( summary( fit3slsd[[ 3 ]]$e1w$eq[[ 2 ]], useDfSys = FALSE ) ),
   digits = 6 ) )

print( round( coef( summary( fit3slsi[[ 4 ]]$e2 ) ), digits = 6 ) )
print( round( coef( summary( fit3slsi[[ 5 ]]$e2$eq[[ 1 ]] ) ), digits = 6 ) )

print( round( coef( summary( fit3slsd[[ 5 ]]$e3e, useDfSys = FALSE ) ),
   digits = 6 ) )
print( round( coef( summary( fit3slsd[[ 5 ]]$e3e, useDfSys = FALSE ),
   modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( summary( fit3slsd[[ 1 ]]$e3e$eq[[ 2 ]], useDfSys = FALSE ) ),
   digits = 6 ) )

print( round( coef( summary( fit3slsi[[ 4 ]]$e3w, useDfSys = FALSE ) ),
   digits = 6 ) )
print( round( coef( summary( fit3slsi[[ 4 ]]$e3w, useDfSys = FALSE ),
   modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( summary( fit3slsi[[ 5 ]]$e3w$eq[[ 2 ]], useDfSys = FALSE ) ),
   digits = 6 ) )

print( round( coef( summary( fit3sls[[ 1 ]]$e4 ) ), digits = 6 ) )
print( round( coef( summary( fit3sls[[ 2 ]]$e4$eq[[ 1 ]] ) ), digits = 6 ) )

print( round( coef( summary( fit3slsi[[ 2 ]]$e5e ) ), digits = 6 ) )
print( round( coef( summary( fit3slsi[[ 2 ]]$e5e ), modified.regMat = TRUE ),
   digits = 6 ) )
print( round( coef( summary( fit3slsi[[ 3 ]]$e5e$eq[[ 2 ]] ) ), digits = 6 ) )

print( round( coef( summary( fit3sls[[ 2 ]]$e5we ) ), digits = 6 ) )
print( round( coef( summary( fit3sls[[ 2 ]]$e5we ), modified.regMat = TRUE ),
   digits = 6 ) )
print( round( coef( summary( fit3sls[[ 3 ]]$e5we$eq[[ 2 ]] ) ), digits = 6 ) )


## *********** variance covariance matrix of the coefficients *******
print( round( vcov( fit3sls[[ 3 ]]$e1c ), digits = 5 ) )
print( round( vcov( fit3sls[[ 4 ]]$e1c$eq[[ 2 ]] ), digits = 5 ) )

print( round( vcov( fit3sls[[ 4 ]]$e2 ), digits = 5 ) )
print( round( vcov( fit3sls[[ 5 ]]$e2$eq[[ 1 ]] ), digits = 5 ) )

print( round( vcov( fit3sls[[ 5 ]]$e3e ), digits = 5 ) )
print( round( vcov( fit3sls[[ 5 ]]$e3e, modified.regMat = TRUE ), digits = 5 ) )
print( round( vcov( fit3sls[[ 1 ]]$e3e$eq[[ 2 ]] ), digits = 5 ) )

print( round( vcov( fit3sls[[ 1 ]]$e4 ), digits = 5 ) )
print( round( vcov( fit3sls[[ 2 ]]$e4$eq[[ 1 ]] ), digits = 5 ) )

print( round( vcov( fit3sls[[ 3 ]]$e4wSym ), digits = 5 ) )
print( round( vcov( fit3sls[[ 4 ]]$e4wSym$eq[[ 1 ]] ), digits = 5 ) )

print( round( vcov( fit3sls[[ 2 ]]$e5e ), digits = 5 ) )
print( round( vcov( fit3sls[[ 2 ]]$e5e, modified.regMat = TRUE ), digits = 5 ) )
print( round( vcov( fit3sls[[ 3 ]]$e5e$eq[[ 2 ]] ), digits = 5 ) )

print( round( vcov( fit3slsi[[ 4 ]]$e1e ), digits = 5 ) )
print( round( vcov( fit3slsi[[ 3 ]]$e1e$eq[[ 1 ]] ), digits = 5 ) )

print( round( vcov( fit3slsi[[ 5 ]]$e1we ), digits = 5 ) )
print( round( vcov( fit3slsi[[ 1 ]]$e1we$eq[[ 2 ]] ), digits = 5 ) )

print( round( vcov( fit3slsi[[ 5 ]]$e2e ), digits = 5 ) )
print( round( vcov( fit3slsi[[ 4 ]]$e2e$eq[[ 2 ]] ), digits = 5 ) )

print( round( vcov( fit3slsi[[ 1 ]]$e3 ), digits = 5 ) )
print( round( vcov( fit3slsi[[ 1 ]]$e3, modified.regMat = TRUE ), digits = 5 ) )
print( round( vcov( fit3slsi[[ 5 ]]$e3$eq[[ 1 ]] ), digits = 5 ) )

print( round( vcov( fit3slsi[[ 2 ]]$e4e ), digits = 5 ) )
print( round( vcov( fit3slsi[[ 1 ]]$e4e$eq[[ 2 ]] ), digits = 5 ) )

print( round( vcov( fit3slsi[[ 3 ]]$e5 ), digits = 5 ) )
print( round( vcov( fit3slsi[[ 3 ]]$e5, modified.regMat = TRUE ), digits = 5 ) )
print( round( vcov( fit3slsi[[ 2 ]]$e5$eq[[ 1 ]] ), digits = 5 ) )

print( round( vcov( fit3slsi[[ 5 ]]$e5w ), digits = 5 ) )
print( round( vcov( fit3slsi[[ 5 ]]$e5w, modified.regMat = TRUE ), digits = 5 ) )
print( round( vcov( fit3slsi[[ 4 ]]$e5w$eq[[ 1 ]] ), digits = 5 ) )

print( round( vcov( fit3slsd[[ 5 ]]$e1c ), digits = 5 ) )
print( round( vcov( fit3slsd[[ 2 ]]$e1c$eq[[ 2 ]] ), digits = 5 ) )

print( round( vcov( fit3slsd[[ 1 ]]$e2 ), digits = 5 ) )
print( round( vcov( fit3slsd[[ 3 ]]$e2$eq[[ 1 ]] ), digits = 5 ) )

print( round( vcov( fit3slsd[[ 5 ]]$e2we ), digits = 5 ) )
print( round( vcov( fit3slsd[[ 3 ]]$e2we$eq[[ 1 ]] ), digits = 5 ) )

print( round( vcov( fit3slsd[[ 2 ]]$e3 ), digits = 5 ) )
print( round( vcov( fit3slsd[[ 2 ]]$e3, modified.regMat = TRUE ), digits = 5 ) )
print( round( vcov( fit3slsd[[ 4 ]]$e3$eq[[ 2 ]] ), digits = 5 ) )

print( round( vcov( fit3slsd[[ 3 ]]$e4 ), digits = 5 ) )
print( round( vcov( fit3slsd[[ 5 ]]$e4$eq[[ 1 ]] ), digits = 5 ) )

print( round( vcov( fit3slsd[[ 4 ]]$e5e ), digits = 5 ) )
print( round( vcov( fit3slsd[[ 4 ]]$e5e, modified.regMat = TRUE ), digits = 5 ) )
print( round( vcov( fit3slsd[[ 1 ]]$e5e$eq[[ 2 ]] ), digits = 5 ) )


## *********** confidence intervals of coefficients *************
print( confint( fit3sls[[ 1 ]]$e1c, useDfSys = TRUE ) )
print( confint( fit3sls[[ 1 ]]$e1c$eq[[ 1 ]], level = 0.9, useDfSys = TRUE ) )

print( confint( fit3sls[[ 2 ]]$e2e, level = 0.9, useDfSys = TRUE ) )
print( confint( fit3sls[[ 2 ]]$e2e$eq[[ 2 ]], level = 0.99, useDfSys = TRUE ) )

print( confint( fit3sls[[ 3 ]]$e3, level = 0.99 ) )
print( confint( fit3sls[[ 3 ]]$e3$eq[[ 1 ]], level = 0.5 ) )

print( confint( fit3sls[[ 5 ]]$e3we, level = 0.99 ) )
print( confint( fit3sls[[ 5 ]]$e3we$eq[[ 1 ]], level = 0.5 ) )

print( confint( fit3sls[[ 4 ]]$e4e, level = 0.5, useDfSys = TRUE ) )
print( confint( fit3sls[[ 4 ]]$e4e$eq[[ 2 ]], level = 0.25, useDfSys = TRUE ) )

print( confint( fit3sls[[ 5 ]]$e5, level = 0.25 ) )
print( confint( fit3sls[[ 5 ]]$e5$eq[[ 1 ]], level = 0.975 ) )

print( confint( fit3slsi[[ 2 ]]$e3e, level = 0.975, useDfSys = TRUE ) )
print( confint( fit3slsi[[ 2 ]]$e3e$eq[[ 1 ]], level = 0.999, useDfSys = TRUE ) )

print( confint( fit3slsi[[ 1 ]]$e5w, level = 0.975, useDfSys = TRUE ) )
print( confint( fit3slsi[[ 1 ]]$e5w$eq[[ 1 ]], level = 0.999, useDfSys = TRUE ) )

print( confint( fit3slsd[[ 3 ]]$e4, level = 0.999 ) )
print( confint( fit3slsd[[ 3 ]]$e4$eq[[ 2 ]] ) )

print( confint( fit3slsd[[ 2 ]]$e4w, level = 0.999 ) )
print( confint( fit3slsd[[ 2 ]]$e4w$eq[[ 2 ]] ) )


## *********** fitted values *************
print( fitted( fit3sls[[ 2 ]]$e1c ) )
print( fitted( fit3sls[[ 2 ]]$e1c$eq[[ 1 ]] ) )

print( fitted( fit3sls[[ 1 ]]$e1wc ) )
print( fitted( fit3sls[[ 1 ]]$e1wc$eq[[ 1 ]] ) )

print( fitted( fit3sls[[ 3 ]]$e2e ) )
print( fitted( fit3sls[[ 3 ]]$e2e$eq[[ 2 ]] ) )

print( fitted( fit3sls[[ 4 ]]$e3 ) )
print( fitted( fit3sls[[ 4 ]]$e3$eq[[ 1 ]] ) )

print( fitted( fit3sls[[ 5 ]]$e4e ) )
print( fitted( fit3sls[[ 5 ]]$e4e$eq[[ 2 ]] ) )

print( fitted( fit3sls[[ 1 ]]$e5 ) )
print( fitted( fit3sls[[ 1 ]]$e5$eq[[ 1 ]] ) )

print( fitted( fit3slsi[[ 3 ]]$e3e ) )
print( fitted( fit3slsi[[ 3 ]]$e3e$eq[[ 1 ]] ) )

print( fitted( fit3slsd[[ 4 ]]$e4 ) )
print( fitted( fit3slsd[[ 4 ]]$e4$eq[[ 2 ]] ) )

print( fitted( fit3slsd[[ 2 ]]$e3w ) )
print( fitted( fit3slsd[[ 2 ]]$e3w$eq[[ 2 ]] ) )


## *********** predicted values *************
predictData <- Kmenta
predictData$consump <- NULL
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fit3sls[[ 2 ]]$e1c, se.fit = TRUE, interval = "prediction",
   useDfSys = TRUE ) )
print( predict( fit3sls[[ 2 ]]$e1c$eq[[ 1 ]], se.fit = TRUE, interval = "prediction",
   useDfSys = TRUE ) )

print( predict( fit3sls[[ 3 ]]$e2e, se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData, useDfSys = TRUE ) )
print( predict( fit3sls[[ 3 ]]$e2e$eq[[ 2 ]], se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData, useDfSys = TRUE ) )

print( predict( fit3sls[[ 5 ]]$e2w, se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData, useDfSys = TRUE ) )
print( predict( fit3sls[[ 5 ]]$e2w$eq[[ 2 ]], se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData, useDfSys = TRUE ) )

print( predict( fit3sls[[ 4 ]]$e3, se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )
print( predict( fit3sls[[ 4 ]]$e3$eq[[ 1 ]], se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )

print( predict( fit3sls[[ 5 ]]$e4e, se.fit = TRUE, interval = "confidence",
   level = 0.25, useDfSys = TRUE ) )
print( predict( fit3sls[[ 5 ]]$e4e$eq[[ 2 ]], se.fit = TRUE, interval = "confidence",
   level = 0.25, useDfSys = TRUE ) )

print( predict( fit3sls[[ 1 ]]$e5, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )
print( predict( fit3sls[[ 1 ]]$e5$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )

print( predict( fit3slsi[[ 3 ]]$e3e, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99, useDfSys = TRUE ) )
print( predict( fit3slsi[[ 3 ]]$e3e$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99, useDfSys = TRUE ) )

print( predict( fit3slsi[[ 1 ]]$e5w, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )
print( predict( fit3slsi[[ 1 ]]$e5w$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )

print( predict( fit3slsd[[ 4 ]]$e4, se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )
print( predict( fit3slsd[[ 4 ]]$e4$eq[[ 2 ]], se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )

print( predict( fit3slsd[[ 2 ]]$e3w, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99, useDfSys = TRUE ) )
print( predict( fit3slsd[[ 2 ]]$e3w$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99, useDfSys = TRUE ) )


# predict just one observation
smallData <- data.frame( price = 130, income = 150, farmPrice = 120,
   trend = 25 )

print( predict( fit3sls[[ 3 ]]$e1c, newdata = smallData ) )
print( predict( fit3sls[[ 3 ]]$e1c$eq[[ 1 ]], newdata = smallData ) )

print( predict( fit3sls[[ 4 ]]$e2e, se.fit = TRUE, level = 0.9,
   newdata = smallData ) )
print( predict( fit3sls[[ 5 ]]$e2e$eq[[ 1 ]], se.pred = TRUE, level = 0.99,
   newdata = smallData ) )

print( predict( fit3sls[[ 1]]$e3, interval = "prediction", level = 0.975,
   newdata = smallData ) )
print( predict( fit3sls[[ 1 ]]$e3$eq[[ 1 ]], interval = "confidence", level = 0.8,
   newdata = smallData ) )

print( predict( fit3sls[[ 4]]$e3we, interval = "prediction", level = 0.975,
   newdata = smallData ) )
print( predict( fit3sls[[ 4 ]]$e3we$eq[[ 1 ]], interval = "confidence", level = 0.8,
   newdata = smallData ) )

print( predict( fit3sls[[ 2 ]]$e4e, se.fit = TRUE, interval = "confidence",
   level = 0.999, newdata = smallData ) )
print( predict( fit3sls[[ 2 ]]$e4e$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.75, newdata = smallData ) )

print( predict( fit3sls[[ 3 ]]$e5, se.fit = TRUE, interval = "prediction",
   newdata = smallData ) )
print( predict( fit3sls[[ 3 ]]$e5$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   newdata = smallData ) )

print( predict( fit3slsi[[ 4 ]]$e3e, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = smallData ) )
print( predict( fit3slsd[[ 5 ]]$e4$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, newdata = smallData ) )

print( predict( fit3slsd[[ 2 ]]$e2we, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = smallData ) )
print( predict( fit3slsi[[ 3 ]]$e4we$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, newdata = smallData ) )


## ************ correlation of predicted values ***************
print( correlation.systemfit( fit3sls[[ 1 ]]$e1c, 2, 1 ) )

print( correlation.systemfit( fit3sls[[ 2 ]]$e2e, 1, 2 ) )

print( correlation.systemfit( fit3sls[[ 5 ]]$e2w, 2, 1 ) )

print( correlation.systemfit( fit3sls[[ 3 ]]$e3, 2, 1 ) )

print( correlation.systemfit( fit3sls[[ 4 ]]$e4e, 1, 2 ) )

print( correlation.systemfit( fit3sls[[ 5 ]]$e5, 2, 1 ) )

print( correlation.systemfit( fit3slsi[[ 2 ]]$e3e, 1, 2 ) )

print( correlation.systemfit( fit3slsi[[ 4 ]]$e5w, 1, 2 ) )

print( correlation.systemfit( fit3slsd[[ 3 ]]$e4, 2, 1 ) )


## ************ Log-Likelihood values ***************
print( logLik( fit3sls[[ 1 ]]$e1c ) )
print( logLik( fit3sls[[ 1 ]]$e1c, residCovDiag = TRUE ) )

print( logLik( fit3sls[[ 2 ]]$e2e ) )
print( logLik( fit3sls[[ 2 ]]$e2e, residCovDiag = TRUE ) )

print( logLik( fit3sls[[ 3 ]]$e3 ) )
print( logLik( fit3sls[[ 3 ]]$e3, residCovDiag = TRUE ) )

print( logLik( fit3sls[[ 4 ]]$e4e ) )
print( logLik( fit3sls[[ 4 ]]$e4e, residCovDiag = TRUE ) )

print( logLik( fit3sls[[ 2 ]]$e4wSym ) )
print( logLik( fit3sls[[ 2 ]]$e4wSym, residCovDiag = TRUE ) )

print( logLik( fit3sls[[ 5 ]]$e5 ) )
print( logLik( fit3sls[[ 5 ]]$e5, residCovDiag = TRUE ) )

print( logLik( fit3slsi[[ 2 ]]$e3e ) )
print( logLik( fit3slsi[[ 2 ]]$e3e, residCovDiag = TRUE ) )

print( logLik( fit3slsi[[ 1 ]]$e1we ) )
print( logLik( fit3slsi[[ 1 ]]$e1we, residCovDiag = TRUE ) )

print( logLik( fit3slsd[[ 3 ]]$e4 ) )
print( logLik( fit3slsd[[ 3 ]]$e4, residCovDiag = TRUE ) )

print( logLik( fit3slsd[[ 5 ]]$e2we ) )
print( logLik( fit3slsd[[ 5 ]]$e2we, residCovDiag = TRUE ) )


## ************** F tests ****************
# testing first restriction
print( linearHypothesis( fit3sls[[ 1 ]]$e1, restrm ) )
linearHypothesis( fit3sls[[ 1 ]]$e1, restrict )

print( linearHypothesis( fit3sls[[ 2 ]]$e1e, restrm ) )
linearHypothesis( fit3sls[[ 2 ]]$e1e, restrict )

print( linearHypothesis( fit3sls[[ 3 ]]$e1c, restrm ) )
linearHypothesis( fit3sls[[ 3 ]]$e1c, restrict )

print( linearHypothesis( fit3slsi[[ 4 ]]$e1, restrm ) )
linearHypothesis( fit3slsi[[ 4 ]]$e1, restrict )

print( linearHypothesis( fit3slsd[[ 5 ]]$e1e, restrm ) )
linearHypothesis( fit3slsd[[ 5 ]]$e1e, restrict )

print( linearHypothesis( fit3slsd[[ 2 ]]$e1w, restrm ) )
linearHypothesis( fit3slsd[[ 2 ]]$e1w, restrict )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
restrictOnly2 <- "- demand_price + supply_price = 0.5"
# first restriction not imposed
print( linearHypothesis( fit3sls[[ 5 ]]$e1c, restrOnly2m, restrOnly2q ) )
linearHypothesis( fit3sls[[ 5 ]]$e1c, restrictOnly2 )

print( linearHypothesis( fit3slsi[[ 1 ]]$e1e, restrOnly2m, restrOnly2q ) )
linearHypothesis( fit3slsi[[ 1 ]]$e1e, restrictOnly2 )

print( linearHypothesis( fit3slsi[[ 3 ]]$e1we, restrOnly2m, restrOnly2q ) )
linearHypothesis( fit3slsi[[ 3 ]]$e1we, restrictOnly2 )

print( linearHypothesis( fit3slsd[[ 2 ]]$e1, restrOnly2m, restrOnly2q ) )
linearHypothesis( fit3slsd[[ 2 ]]$e1, restrictOnly2 )

# first restriction imposed
print( linearHypothesis( fit3sls[[ 4 ]]$e2, restrOnly2m, restrOnly2q ) )
linearHypothesis( fit3sls[[ 4 ]]$e2, restrictOnly2 )

print( linearHypothesis( fit3sls[[ 4 ]]$e3, restrOnly2m, restrOnly2q ) )
linearHypothesis( fit3sls[[ 4 ]]$e3, restrictOnly2 )

print( linearHypothesis( fit3sls[[ 1 ]]$e2w, restrOnly2m, restrOnly2q ) )
linearHypothesis( fit3sls[[ 1 ]]$e2w, restrictOnly2 )

print( linearHypothesis( fit3sls[[ 1 ]]$e3we, restrOnly2m, restrOnly2q ) )
linearHypothesis( fit3sls[[ 1 ]]$e3we, restrictOnly2 )

print( linearHypothesis( fit3slsi[[ 5 ]]$e2e, restrOnly2m, restrOnly2q ) )
linearHypothesis( fit3slsi[[ 5 ]]$e2e, restrictOnly2 )

print( linearHypothesis( fit3slsi[[ 5 ]]$e3e, restrOnly2m, restrOnly2q ) )
linearHypothesis( fit3slsi[[ 5 ]]$e3e, restrictOnly2 )

print( linearHypothesis( fit3slsd[[ 1 ]]$e2, restrOnly2m, restrOnly2q ) )
linearHypothesis( fit3slsd[[ 1 ]]$e2, restrictOnly2 )

print( linearHypothesis( fit3slsd[[ 1 ]]$e3, restrOnly2m, restrOnly2q ) )
linearHypothesis( fit3slsd[[ 1 ]]$e3, restrictOnly2 )

# testing both of the restrictions
print( linearHypothesis( fit3sls[[ 2 ]]$e1e, restr2m, restr2q ) )
linearHypothesis( fit3sls[[ 2 ]]$e1e, restrict2 )

print( linearHypothesis( fit3slsi[[ 3 ]]$e1, restr2m, restr2q ) )
linearHypothesis( fit3slsi[[ 3 ]]$e1, restrict2 )

print( linearHypothesis( fit3slsd[[ 4 ]]$e1e, restr2m, restr2q ) )
linearHypothesis( fit3slsd[[ 4 ]]$e1e, restrict2 )

print( linearHypothesis( fit3slsd[[ 5 ]]$e1w, restr2m, restr2q ) )
linearHypothesis( fit3slsd[[ 5 ]]$e1w, restrict2 )


## ************** Wald tests ****************
# testing first restriction
print( linearHypothesis( fit3sls[[ 1 ]]$e1, restrm, test = "Chisq" ) )
linearHypothesis( fit3sls[[ 1 ]]$e1, restrict, test = "Chisq" )

print( linearHypothesis( fit3sls[[ 2 ]]$e1e, restrm, test = "Chisq" ) )
linearHypothesis( fit3sls[[ 2 ]]$e1e, restrict, test = "Chisq" )

print( linearHypothesis( fit3sls[[ 3 ]]$e1c, restrm, test = "Chisq" ) )
linearHypothesis( fit3sls[[ 3 ]]$e1c, restrict, test = "Chisq" )

print( linearHypothesis( fit3slsi[[ 4 ]]$e1, restrm, test = "Chisq" ) )
linearHypothesis( fit3slsi[[ 4 ]]$e1, restrict, test = "Chisq" )

print( linearHypothesis( fit3slsi[[ 2 ]]$e1we, restrm, test = "Chisq" ) )
linearHypothesis( fit3slsi[[ 2 ]]$e1we, restrict, test = "Chisq" )

print( linearHypothesis( fit3slsd[[ 5 ]]$e1e, restrm, test = "Chisq" ) )
linearHypothesis( fit3slsd[[ 5 ]]$e1e, restrict, test = "Chisq" )

# testing second restriction
# first restriction not imposed
print( linearHypothesis( fit3sls[[ 5 ]]$e1c, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fit3sls[[ 5 ]]$e1c, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fit3sls[[ 3 ]]$e1wc, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fit3sls[[ 3 ]]$e1wc, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fit3slsi[[ 1 ]]$e1e, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fit3slsi[[ 1 ]]$e1e, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fit3slsd[[ 2 ]]$e1, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fit3slsd[[ 2 ]]$e1, restrictOnly2, test = "Chisq" )

# first restriction imposed
print( linearHypothesis( fit3sls[[ 4 ]]$e2, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fit3sls[[ 4 ]]$e2, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fit3sls[[ 4 ]]$e3, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fit3sls[[ 4 ]]$e3, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fit3slsi[[ 5 ]]$e2e, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fit3slsi[[ 5 ]]$e2e, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fit3slsi[[ 5 ]]$e3e, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fit3slsi[[ 5 ]]$e3e, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fit3slsd[[ 1 ]]$e2, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fit3slsd[[ 1 ]]$e2, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fit3slsd[[ 1 ]]$e3, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fit3slsd[[ 1 ]]$e3, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fit3slsd[[ 2 ]]$e2we, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fit3slsd[[ 2 ]]$e2we, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fit3slsd[[ 3 ]]$e3w, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fit3slsd[[ 3 ]]$e3w, restrictOnly2, test = "Chisq" )

# testing both of the restrictions
print( linearHypothesis( fit3sls[[ 2 ]]$e1e, restr2m, restr2q, test = "Chisq" ) )
linearHypothesis( fit3sls[[ 2 ]]$e1e, restrict2, test = "Chisq" )

print( linearHypothesis( fit3sls[[ 5 ]]$e1wc, restr2m, restr2q, test = "Chisq" ) )
linearHypothesis( fit3sls[[ 5 ]]$e1wc, restrict2, test = "Chisq" )

print( linearHypothesis( fit3slsi[[ 3 ]]$e1, restr2m, restr2q, test = "Chisq" ) )
linearHypothesis( fit3slsi[[ 3 ]]$e1, restrict2, test = "Chisq" )

print( linearHypothesis( fit3slsd[[ 4 ]]$e1e, restr2m, restr2q, test = "Chisq" ) )
linearHypothesis( fit3slsd[[ 4 ]]$e1e, restrict2, test = "Chisq" )


## *********** model frame *************
print( mf <- model.frame( fit3sls[[ 3 ]]$e1c ) )
print( mf1 <- model.frame( fit3sls[[ 3 ]]$e1c$eq[[ 1 ]] ) )
print( attributes( mf1 )$terms )
print( mf2 <- model.frame( fit3sls[[ 3 ]]$e1c$eq[[ 2 ]] ) )
print( attributes( mf2 )$terms )

print( all.equal( mf, model.frame( fit3sls[[ 3 ]]$e1wc ) ) )
print( all.equal( mf2, model.frame( fit3sls[[ 3 ]]$e1wc$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fit3sls[[ 4 ]]$e2e ) ) )
print( all.equal( mf2, model.frame( fit3sls[[ 4 ]]$e2e$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fit3sls[[ 5 ]]$e3 ) ) )
print( all.equal( mf1, model.frame( fit3sls[[ 5 ]]$e3$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fit3sls[[ 1 ]]$e4e ) ) )
print( all.equal( mf2, model.frame( fit3sls[[ 1 ]]$e4e$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fit3sls[[ 2 ]]$e5 ) ) )
print( all.equal( mf1, model.frame( fit3sls[[ 3 ]]$e5$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fit3slsi[[ 4 ]]$e3e ) ) )
print( all.equal( mf1, model.frame( fit3slsi[[ 4 ]]$e3e$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fit3slsd[[ 5 ]]$e4 ) ) )
print( all.equal( mf2, model.frame( fit3slsd[[ 5 ]]$e4$eq[[ 2 ]] ) ) )

fit3sls[[ 3 ]]$e1c$eq[[ 1 ]]$modelInst
fit3sls[[ 3 ]]$e1c$eq[[ 2 ]]$modelInst

fit3sls[[ 1 ]]$e3$eq[[ 1 ]]$modelInst
fit3sls[[ 1 ]]$e3$eq[[ 2 ]]$modelInst

fit3slsd[[ 5 ]]$e4$eq[[ 1 ]]$modelInst
fit3slsd[[ 5 ]]$e4$eq[[ 2 ]]$modelInst


## **************** model matrix ************************
# with x (returnModelMatrix) = TRUE
print( !is.null( fit3sls[[ 4 ]]$e1c$eq[[ 1 ]]$x ) )
print( mm <- model.matrix( fit3sls[[ 4 ]]$e1c ) )
print( mm1 <- model.matrix( fit3sls[[ 4 ]]$e1c$eq[[ 1 ]] ) )
print( mm2 <- model.matrix( fit3sls[[ 4 ]]$e1c$eq[[ 2 ]] ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fit3sls[[ 4 ]]$e1wc ) ) )
print( all.equal( mm1, model.matrix( fit3sls[[ 4 ]]$e1wc$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit3sls[[ 4 ]]$e1wc$eq[[ 2 ]] ) ) )
print( !is.null( fit3sls[[ 4 ]]$e1wc$eq[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fit3sls[[ 5 ]]$e2$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fit3sls[[ 5 ]]$e2 ) ) )
print( all.equal( mm1, model.matrix( fit3sls[[ 5 ]]$e2$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit3sls[[ 5 ]]$e2$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fit3sls[[ 5 ]]$e2e ) ) )
print( all.equal( mm1, model.matrix( fit3sls[[ 5 ]]$e2e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit3sls[[ 5 ]]$e2e$eq[[ 2 ]] ) ) )
print( !is.null( fit3sls[[ 5 ]]$e1wc$e2e[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fit3sls[[ 1 ]]$e3e$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fit3sls[[ 1 ]]$e3e ) ) )
print( all.equal( mm1, model.matrix( fit3sls[[ 1 ]]$e3e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit3sls[[ 1 ]]$e3e$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fit3sls[[ 1 ]]$e3 ) ) )
print( all.equal( mm1, model.matrix( fit3sls[[ 1 ]]$e3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit3sls[[ 1 ]]$e3$eq[[ 2 ]] ) ) )
print( !is.null( fit3sls[[ 1 ]]$e3$eq[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fit3slsi[[ 2 ]]$e4$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fit3slsi[[ 2 ]]$e4 ) ) )
print( all.equal( mm1, model.matrix( fit3slsi[[ 2 ]]$e4$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit3slsi[[ 2 ]]$e4$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fit3slsi[[ 2 ]]$e4we ) ) )
print( all.equal( mm1, model.matrix( fit3slsi[[ 2 ]]$e4we$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit3slsi[[ 2 ]]$e4we$eq[[ 2 ]] ) ) )
print( !is.null( fit3slsi[[ 2 ]]$e1wc$e4we[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fit3slsi[[ 5 ]]$e5w$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fit3slsi[[ 5 ]]$e5w ) ) )
print( all.equal( mm1, model.matrix( fit3slsi[[ 5 ]]$e5w$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit3slsi[[ 5 ]]$e5w$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fit3slsi[[ 5 ]]$e5 ) ) )
print( all.equal( mm1, model.matrix( fit3slsi[[ 5 ]]$e5$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit3slsi[[ 5 ]]$e5$eq[[ 2 ]] ) ) )
print( !is.null( fit3slsi[[ 5 ]]$e5$eq[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fit3slsd[[ 3 ]]$e5e$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fit3slsd[[ 3 ]]$e5e ) ) )
print( all.equal( mm1, model.matrix( fit3slsd[[ 3 ]]$e5e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit3slsd[[ 3 ]]$e5e$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fit3slsd[[ 3 ]]$e5we ) ) )
print( all.equal( mm1, model.matrix( fit3slsd[[ 3 ]]$e5we$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit3slsd[[ 3 ]]$e5we$eq[[ 2 ]] ) ) )
print( !is.null( fit3sls[[ 3 ]]$e5we$eq[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fit3slsd[[ 2 ]]$e3w$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fit3slsd[[ 2 ]]$e3w ) ) )
print( all.equal( mm1, model.matrix( fit3slsd[[ 2 ]]$e3w$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit3slsd[[ 2 ]]$e3w$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fit3slsd[[ 2 ]]$e3 ) ) )
print( all.equal( mm1, model.matrix( fit3slsd[[ 2 ]]$e3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit3slsd[[ 2 ]]$e3$eq[[ 2 ]] ) ) )
print( !is.null( fit3slsd[[ 2 ]]$e3$eq[[ 1 ]]$x ) )

# matrices of instrumental variables
model.matrix( fit3sls[[ 1 ]]$e1c, which = "z" )
model.matrix( fit3sls[[ 3 ]]$e1c$eq[[ 1 ]], which = "z" )
model.matrix( fit3sls[[ 4 ]]$e1c$eq[[ 2 ]], which = "z" )

# matrices of fitted regressors
model.matrix( fit3slsd[[ 1 ]]$e3w, which = "xHat" )
model.matrix( fit3slsd[[ 3 ]]$e3w$eq[[ 1 ]], which = "xHat" )
model.matrix( fit3slsd[[ 4 ]]$e3w$eq[[ 2 ]], which = "xHat" )


## **************** formulas ************************
formula( fit3sls[[ 2 ]]$e1c )
formula( fit3sls[[ 2 ]]$e1c$eq[[ 1 ]] )

formula( fit3sls[[ 3 ]]$e2e )
formula( fit3sls[[ 3 ]]$e2e$eq[[ 2 ]] )

formula( fit3sls[[ 4 ]]$e3 )
formula( fit3sls[[ 4 ]]$e3$eq[[ 1 ]] )

formula( fit3sls[[ 5 ]]$e4e )
formula( fit3sls[[ 5 ]]$e4e$eq[[ 2 ]] )

formula( fit3sls[[ 1 ]]$e5 )
formula( fit3sls[[ 1 ]]$e5$eq[[ 1 ]] )

formula( fit3slsi[[ 3 ]]$e3e )
formula( fit3slsi[[ 3 ]]$e3e$eq[[ 1 ]] )

formula( fit3slsd[[ 4 ]]$e4 )
formula( fit3slsd[[ 4 ]]$e4$eq[[ 2 ]] )

formula( fit3slsd[[ 2 ]]$e1w )
formula( fit3slsd[[ 2 ]]$e1w$eq[[ 1 ]] )


## **************** model terms *******************
terms( fit3sls[[ 2 ]]$e1c )
terms( fit3sls[[ 2 ]]$e1c$eq[[ 1 ]] )

terms( fit3sls[[ 3 ]]$e2e )
terms( fit3sls[[ 3 ]]$e2e$eq[[ 2 ]] )

terms( fit3sls[[ 4 ]]$e3 )
terms( fit3sls[[ 4 ]]$e3$eq[[ 1 ]] )

terms( fit3sls[[ 5 ]]$e4e )
terms( fit3sls[[ 5 ]]$e4e$eq[[ 2 ]] )

terms( fit3sls[[ 1 ]]$e5 )
terms( fit3sls[[ 1 ]]$e5$eq[[ 1 ]] )

terms( fit3sls[[ 2 ]]$e4wSym )
terms( fit3sls[[ 2 ]]$e4wSym$eq[[ 1 ]] )

terms( fit3slsi[[ 3 ]]$e3e )
terms( fit3slsi[[ 3 ]]$e3e$eq[[ 1 ]] )

terms( fit3slsd[[ 4 ]]$e4 )
terms( fit3slsd[[ 4 ]]$e4$eq[[ 2 ]] )

terms( fit3slsd[[ 5 ]]$e5we )
terms( fit3slsd[[ 5 ]]$e5we$eq[[ 2 ]] )


## **************** terms of instruments *******************
fit3sls[[ 2 ]]$e1c$eq[[ 1 ]]$termsInst

fit3sls[[ 3 ]]$e2e$eq[[ 2 ]]$termsInst

fit3sls[[ 4 ]]$e3$eq[[ 1 ]]$termsInst

fit3sls[[ 5 ]]$e4e$eq[[ 2 ]]$termsInst

fit3sls[[ 1 ]]$e5$eq[[ 1 ]]$termsInst

fit3sls[[ 2 ]]$e4wSym$eq[[ 1 ]]$termsInst

fit3slsi[[ 3 ]]$e3e$eq[[ 1 ]]$termsInst

fit3slsd[[ 4 ]]$e4$eq[[ 2 ]]$termsInst

fit3slsd[[ 5 ]]$e5we$eq[[ 2 ]]$termsInst


## **************** estfun ************************
library( "sandwich" )

estfun( fit3sls[[ 1 ]]$e1 )
round( colSums( estfun( fit3sls[[ 1 ]]$e1 ) ), digits = 7 )

estfun( fit3sls[[ 2 ]]$e1e )
round( colSums( estfun( fit3sls[[ 2 ]]$e1e ) ), digits = 7 )

estfun( fit3sls[[ 3 ]]$e1c )
round( colSums( estfun( fit3sls[[ 3 ]]$e1c ) ), digits = 7 )

estfun( fit3sls[[ 4 ]]$e1wc )

round( colSums( estfun( fit3sls[[ 5 ]]$e1wc ) ), digits = 7 )
round( colSums( estfun( fit3sls[[ 5 ]]$e1wc, residFit = FALSE ) ), digits = 7 )

round( colSums( estfun( fit3sls[[ 4 ]]$e1wc ) ), digits = 7 )
round( colSums( estfun( fit3sls[[ 4 ]]$e1wc, residFit = FALSE ) ), digits = 7 )

round( colSums( estfun( fit3sls[[ 3 ]]$e1wc ) ), digits = 7 )
round( colSums( estfun( fit3sls[[ 3 ]]$e1wc, residFit = FALSE ) ), digits = 7 )

round( colSums( estfun( fit3sls[[ 2 ]]$e1wc ) ), digits = 7 )
round( colSums( estfun( fit3sls[[ 2 ]]$e1wc, residFit = FALSE ) ), digits = 7 )

round( colSums( estfun( fit3sls[[ 1 ]]$e1wc ) ), digits = 7 )
round( colSums( estfun( fit3sls[[ 1 ]]$e1wc, residFit = FALSE ) ), digits = 7 )

estfun( fit3slsd[[ 5 ]]$e1w )
estfun( fit3slsd[[ 5 ]]$e1w, residFit = FALSE )

round( colSums( estfun( fit3slsd[[ 5 ]]$e1w ) ), digits = 7 )
round( colSums( estfun( fit3slsd[[ 5 ]]$e1w, residFit = FALSE ) ), digits = 7 )

round( colSums( estfun( fit3slsd[[ 4 ]]$e1w ) ), digits = 7 )
round( colSums( estfun( fit3slsd[[ 4 ]]$e1w, residFit = FALSE ) ), digits = 7 )

round( colSums( estfun( fit3slsd[[ 3 ]]$e1w ) ), digits = 7 )
round( colSums( estfun( fit3slsd[[ 3 ]]$e1w, residFit = FALSE ) ), digits = 7 )

round( colSums( estfun( fit3slsd[[ 2 ]]$e1w ) ), digits = 7 )
round( colSums( estfun( fit3slsd[[ 2 ]]$e1w, residFit = FALSE ) ), digits = 7 )

round( colSums( estfun( fit3slsd[[ 1 ]]$e1w ) ), digits = 7 )
round( colSums( estfun( fit3slsd[[ 1 ]]$e1w, residFit = FALSE ) ), digits = 7 )


## **************** bread ************************
bread( fit3sls[[ 1 ]]$e1 )

bread( fit3sls[[ 2 ]]$e1e )

bread( fit3sls[[ 3 ]]$e1c )

bread( fit3sls[[ 4 ]]$e1wc )

bread( fit3slsd[[ 5 ]]$e1w )

bread( fit3slsd[[ 4 ]]$e1w )

bread( fit3slsd[[ 3 ]]$e1w )

bread( fit3slsd[[ 2 ]]$e1w )

bread( fit3slsd[[ 1 ]]$e1w )
