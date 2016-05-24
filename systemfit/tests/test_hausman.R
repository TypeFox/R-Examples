library( "systemfit" )
options( digits = 5 )

data( "Kmenta" )
useMatrix <- FALSE

eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
inst <- ~ income + farmPrice + trend
eqSystem <- list( demand = eqDemand, supply = eqSupply )
restrm <- matrix(0,1,7)  # restriction matrix "R"
restrm[1,3] <-  1
restrm[1,7] <- -1
restr2m <- matrix(0,2,7)  # restriction matrix "R" 2
restr2m[1,3] <-  1
restr2m[1,7] <- -1
restr2m[2,2] <- -1
restr2m[2,5] <-  1
restr2q <- c( 0, 0.5 )  # restriction vector "q" 2
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


## ******************* unrestricted estimation *****************
## ******************** default estimation *********************
fit2sls1 <- systemfit( eqSystem, "2SLS", data = Kmenta, inst = inst,
   useMatrix = useMatrix )
fit3sls1 <- systemfit( eqSystem, "3SLS", data = Kmenta, inst = inst,
   useMatrix = useMatrix )
print( hausman.systemfit( fit2sls1, fit3sls1 ) )

## ************** 2SLS estimation with singleEqSigma = FALSE *****************
fit2sls1s <- systemfit( eqSystem, "2SLS", data = Kmenta, inst = inst,
   singleEqSigma = FALSE, useMatrix = useMatrix )
print( hausman.systemfit( fit2sls1s, fit3sls1 ) )

## ******************* estimations with methodResidCov = 0 *****************
fit2sls1r <- systemfit( eqSystem, "2SLS", data = Kmenta, inst = inst,
   methodResidCov = "noDfCor", useMatrix = useMatrix )
fit3sls1r <- systemfit( eqSystem, "3SLS", data = Kmenta, inst = inst,
   methodResidCov = "noDfCor", useMatrix = useMatrix )
print( hausman.systemfit( fit2sls1r, fit3sls1r ) )


## ********************* estimation with restriction ********************
## *********************** default estimation ***********************
fit2sls2 <- systemfit( eqSystem, "2SLS", data = Kmenta, restrict.matrix = restrm,
   inst = inst, useMatrix = useMatrix )
fit3sls2 <- systemfit( eqSystem, "3SLS", data = Kmenta, restrict.matrix = restrm,
   inst = inst, useMatrix = useMatrix )
# print( hausman.systemfit( fit2sls2, fit3sls2 ) )

## ************* 2SLS estimation with singleEqSigma = TRUE *****************
fit2sls2s <- systemfit( eqSystem, "2SLS", data = Kmenta, restrict.matrix = restrm,
   inst = inst, singleEqSigma = TRUE, useMatrix = useMatrix )
# print( hausman.systemfit( fit2sls2s, fit3sls2 ) )

## ********************* estimations with methodResidCov = 0 **************
fit2sls2r <- systemfit( eqSystem, "2SLS", data = Kmenta, restrict.matrix = restrm,
   inst = inst, methodResidCov = "noDfCor", useMatrix = useMatrix )
fit3sls2r <- systemfit( eqSystem, "3SLS", data = Kmenta, restrict.matrix = restrm,
   inst = inst, methodResidCov = "noDfCor", useMatrix = useMatrix )
# print( hausman.systemfit( fit2sls2r, fit3sls2r ) )


## ****************** estimation with restriction via restrict.regMat ******************
## ********************** default estimation ********************
fit2sls3 <- systemfit( eqSystem, "2SLS", data = Kmenta, restrict.regMat = tc,
   inst = inst, useMatrix = useMatrix )
fit3sls3 <- systemfit( eqSystem, "3SLS", data = Kmenta, restrict.regMat = tc,
   inst = inst, useMatrix = useMatrix )
print( hausman.systemfit( fit2sls3, fit3sls3 ) )

## ******************* estimations with methodResidCov = 0 *******
fit2sls3r <- systemfit( eqSystem, "2SLS", data = Kmenta, restrict.regMat = tc,
   inst = inst, methodResidCov = "noDfCor", useMatrix = useMatrix )
fit3sls3r <- systemfit( eqSystem, "3SLS", data = Kmenta, restrict.regMat = tc,
   inst = inst, methodResidCov = "noDfCor", useMatrix = useMatrix )
print( hausman.systemfit( fit2sls3r, fit3sls3r ) )


## ***************** estimations with 2 restrictions *******************
## *********************** default estimations **************
fit2sls4 <- systemfit( eqSystem, "2SLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, inst = inst, useMatrix = useMatrix )
fit3sls4 <- systemfit( eqSystem, "3SLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, inst = inst, useMatrix = useMatrix )
# print( hausman.systemfit( fit2sls4, fit3sls4 ) )

## ***************** estimations with methodResidCov = 0 **************
fit2sls4r <- systemfit( eqSystem, "2SLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, inst = inst, methodResidCov = "noDfCor",
   useMatrix = useMatrix )
fit3sls4r <- systemfit( eqSystem, "3SLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, inst = inst, methodResidCov = "noDfCor",
   useMatrix = useMatrix )
# print( hausman.systemfit( fit2sls4r, fit3sls4r ) )


## *********** estimations with 2 restrictions via R and restrict.regMat ***************
## ***************** default estimations *******************
fit2sls5 <- systemfit( eqSystem, "2SLS", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, inst = inst,
   useMatrix = useMatrix )
fit3sls5 <- systemfit( eqSystem, "3SLS", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, inst = inst,
   useMatrix = useMatrix )
# print( hausman.systemfit( fit2sls5, fit3sls5 ) )

## ************* estimations with methodResidCov = 0 *********
fit2sls5r <- systemfit( eqSystem, "2SLS", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, inst = inst,
   methodResidCov = "noDfCor", useMatrix = useMatrix )
fit3sls5r <- systemfit( eqSystem, "3SLS", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, inst = inst,
   methodResidCov = "noDfCor", useMatrix = useMatrix )
# print( hausman.systemfit( fit2sls5r, fit3sls5r ) )
