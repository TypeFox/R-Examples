library( systemfit )
options( digits = 3 )

data( "Kmenta" )
useMatrix <- FALSE

demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
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
restrict2i <- c( "demand_income - supply_trend = 0",
   "- demand_price + supply_income = 0.5" )
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

# the standard equations do not converge and lead to a singular weighting matrix
# both in R and in EViews, since both equations have the same endogenous variable
supply2 <- price ~ income + farmPrice + trend
system2 <- list( demand = demand, supply = supply2 )


## *************** SUR estimation ************************
fitsur1 <- systemfit( system, "SUR", data = Kmenta, useMatrix = useMatrix )
print( summary( fitsur1 ) )
nobs( fitsur1 )

## ********************* SUR (EViews-like) *****************
fitsur1e <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "noDfCor",
   useMatrix = useMatrix )
print( summary( fitsur1e, useDfSys = TRUE ) )
nobs( fitsur1e )

## ********************* SUR (methodResidCov="Theil") *****************
fitsur1r2 <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "Theil",
   useMatrix = useMatrix )
print( summary( fitsur1r2 ) )

## *************** SUR (methodResidCov="Theil", useDfSys = TRUE ) ***************
fitsur1e2 <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "Theil",
   x = TRUE, useMatrix = useMatrix )
print( summary( fitsur1e2, useDfSys = TRUE ) )

## ********************* SUR (methodResidCov="max") *****************
fitsur1r3 <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "max",
   useMatrix = useMatrix )
print( summary( fitsur1r3 ) )

## *************** WSUR estimation ************************
fitsur1w <- systemfit( system, "SUR", data = Kmenta, residCovWeighted = TRUE,
   useMatrix = useMatrix )
summary( fitsur1w )
nobs( fitsur1w )

## *************** WSUR (methodResidCov="Theil", useDfSys = TRUE ) ***************
fitsur1we2 <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "Theil",
   residCovWeighted = TRUE, useMatrix = useMatrix )
summary( fitsur1we2, useDfSys = TRUE )


## *************** SUR with cross-equation restriction **************
fitsur2 <- systemfit( system, "SUR", data = Kmenta, restrict.matrix = restrm,
   useMatrix = useMatrix )
print( summary( fitsur2 ) )
nobs( fitsur2 )
# the same with symbolically specified restrictions
fitsur2Sym <- systemfit( system, "SUR", data = Kmenta,
   restrict.matrix = restrict, useMatrix = useMatrix )
all.equal( fitsur2, fitsur2Sym )
nobs( fitsur2Sym )

## *************** SUR with cross-equation restriction (EViews-like) **
fitsur2e <- systemfit( system, "SUR", data = Kmenta, restrict.matrix = restrm,
   methodResidCov = "noDfCor", x = TRUE,
   useMatrix = useMatrix )
print( summary( fitsur2e ) )

## *************** WSUR with cross-equation restriction (EViews-like) **
fitsur2we <- systemfit( system, "SUR", data = Kmenta, restrict.matrix = restrm,
   methodResidCov = "noDfCor", residCovWeighted = TRUE,
   x = TRUE, useMatrix = useMatrix )
summary( fitsur2we )


## *************** SUR with restriction via restrict.regMat *******************
fitsur3 <- systemfit( system, "SUR", data = Kmenta, restrict.regMat = tc,
   useMatrix = useMatrix )
print( summary( fitsur3 ) )
nobs( fitsur3 )

## *************** SUR with restriction via restrict.regMat (EViews-like) **************
fitsur3e <- systemfit( system, "SUR", data = Kmenta, restrict.regMat = tc,
   methodResidCov = "noDfCor", x = TRUE,
   useMatrix = useMatrix )
print( summary( fitsur3e ) )

## *************** WSUR with restriction via restrict.regMat *******************
fitsur3w <- systemfit( system, "SUR", data = Kmenta, restrict.regMat = tc,
   residCovWeighted = TRUE, x = TRUE, useMatrix = useMatrix )
summary( fitsur3w )


## *************** SUR with 2 restrictions ***************************
fitsur4 <- systemfit( system, "SUR", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, useMatrix = useMatrix )
print( summary( fitsur4 ) )
nobs( fitsur4 )
# the same with symbolically specified restrictions
fitsur4Sym <- systemfit( system, "SUR", data = Kmenta,
   restrict.matrix = restrict2, useMatrix = useMatrix )
all.equal( fitsur4, fitsur4Sym )

## *************** SUR with 2 restrictions (EViews-like) **************
fitsur4e <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "noDfCor",
   restrict.matrix = restr2m, restrict.rhs = restr2q, useMatrix = useMatrix )
print( summary( fitsur4e ) )

## *************** SUR with 2 restrictions (methodResidCov = "Theil") **************
fitsur4r2 <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "Theil",
   restrict.matrix = restr2m, restrict.rhs = restr2q, useMatrix = useMatrix )
print( summary( fitsur4r2 ) )

## *************** SUR with 2 restrictions (methodResidCov = "max") **************
fitsur4r3 <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "max",
   restrict.matrix = restr2m, restrict.rhs = restr2q,
   x = TRUE, useMatrix = useMatrix )
print( summary( fitsur4r3 ) )

## *************** WSUR with 2 restrictions (EViews-like) **************
fitsur4we <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "noDfCor",
   restrict.matrix = restr2m, restrict.rhs = restr2q, residCovWeighted = TRUE,
   useMatrix = useMatrix )
summary( fitsur4we )


## *************** SUR with 2 restrictions via R and restrict.regMat ****************
fitsur5 <- systemfit( system, "SUR", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc,
   x = TRUE, useMatrix = useMatrix )
print( summary( fitsur5 ) )
nobs( fitsur5 )
# the same with symbolically specified restrictions
fitsur5Sym <- systemfit( system, "SUR", data = Kmenta,
   restrict.matrix = restrict3, restrict.regMat = tc,
   x = TRUE, useMatrix = useMatrix )
all.equal( fitsur5, fitsur5Sym )

## *************** SUR with 2 restrictions via R and restrict.regMat (EViews-like) **************
fitsur5e <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "noDfCor",
   restrict.matrix = restr3m, restrict.rhs = restr3q, restrict.regMat = tc,
   useMatrix = useMatrix )
print( summary( fitsur5e ) )

## ************ WSUR with 2 restrictions via R and restrict.regMat ************
fitsur5w <- systemfit( system, "SUR", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, residCovWeighted = TRUE,
   useMatrix = useMatrix )
summary( fitsur5w )


## ************** iterated SUR ****************************
fitsuri1 <- systemfit( system2, "SUR", data = Kmenta, maxit = 100,
   useMatrix = useMatrix )
print( summary( fitsuri1 ) )
nobs( fitsuri1 )

## ************** iterated SUR (EViews-like) *****************
fitsuri1e <- systemfit( system2, "SUR", data = Kmenta, methodResidCov = "noDfCor",
   maxit = 100, useMatrix = useMatrix )
print( summary( fitsuri1e, useDfSys = TRUE ) )

## ************** iterated SUR (methodResidCov = "Theil") ****************************
fitsuri1r2 <- systemfit( system2, "SUR", data = Kmenta, maxit = 100,
   methodResidCov = "Theil", useMatrix = useMatrix )
print( summary( fitsuri1r2 ) )

## ************** iterated SUR (methodResidCov="Theil", useDfSys=TRUE) *****************
fitsuri1e2 <- systemfit( system2, "SUR", data = Kmenta, methodResidCov = "Theil",
   maxit = 100, x = TRUE, useMatrix = useMatrix )
print( summary( fitsuri1e2, useDfSys = TRUE ) )

## ************** iterated SUR (methodResidCov = "max") ****************************
fitsuri1r3 <- systemfit( system2, "SUR", data = Kmenta, maxit = 100,
   methodResidCov = "max", useMatrix = useMatrix )
print( summary( fitsuri1r3 ) )

## ************** iterated WSUR (methodResidCov = "max") ****************************
fitsuri1wr3 <- systemfit( system2, "SUR", data = Kmenta, maxit = 100,
   methodResidCov = "max", residCovWeighted = TRUE, useMatrix = useMatrix )
summary( fitsuri1wr3 )


## *********** iterated SUR with restriction *******************
fitsuri2 <- systemfit( system2, "SUR", data = Kmenta, restrict.matrix = restrm,
   maxit = 100, useMatrix = useMatrix )
print( summary( fitsuri2 ) )

## *********** iterated SUR with restriction (EViews-like) ***************
fitsuri2e <- systemfit( system2, "SUR", data = Kmenta, restrict.matrix = restrm,
   methodResidCov = "noDfCor", maxit = 100, x = TRUE,
   useMatrix = useMatrix )
print( summary( fitsuri2e ) )

## *********** iterated WSUR with restriction *******************
fitsuri2w <- systemfit( system2, "SUR", data = Kmenta, restrict.matrix = restrm,
   maxit = 100, residCovWeighted = TRUE, useMatrix = useMatrix )
summary( fitsuri2w )


## *********** iterated SUR with restriction via restrict.regMat ********************
fitsuri3 <- systemfit( system2, "SUR", data = Kmenta, restrict.regMat = tc,
   maxit = 100, useMatrix = useMatrix )
print( summary( fitsuri3 ) )

## *********** iterated SUR with restriction via restrict.regMat (EViews-like) ***************
fitsuri3e <- systemfit( system2, "SUR", data = Kmenta, restrict.regMat = tc,
   methodResidCov = "noDfCor", maxit = 100, x = TRUE,
   useMatrix = useMatrix )
print( summary( fitsuri3e ) )

## *********** iterated WSUR with restriction via restrict.regMat (EViews-like) ***************
fitsuri3we <- systemfit( system2, "SUR", data = Kmenta, restrict.regMat = tc,
   methodResidCov = "noDfCor", maxit = 100, residCovWeighted = TRUE,
   useMatrix = useMatrix )
summary( fitsuri3we )


## *************** iterated SUR with 2 restrictions ***************************
fitsurio4 <- systemfit( system, "SUR", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, maxit = 100, useMatrix = useMatrix )
print( summary( fitsurio4 ) )
fitsuri4 <- systemfit( system2, "SUR", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, maxit = 100, useMatrix = useMatrix )
print( summary( fitsuri4 ) )

## *************** iterated SUR with 2 restrictions (EViews-like) **************
fitsurio4e <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "noDfCor",
   restrict.matrix = restr2m, restrict.rhs = restr2q, maxit = 100,
   useMatrix = useMatrix )
print( summary( fitsurio4e ) )
fitsuri4e <- systemfit( system2, "SUR", data = Kmenta, methodResidCov = "noDfCor",
   restrict.matrix = restr2m, restrict.rhs = restr2q, maxit = 100,
   useMatrix = useMatrix )
print( summary( fitsuri4e ) )

## *************** iterated WSUR with 2 restrictions ***************************
fitsurio4w <- systemfit( system, "SUR", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, maxit = 100, residCovWeighted = TRUE,
   useMatrix = useMatrix )
summary( fitsurio4w )
fitsuri4w <- systemfit( system2, "SUR", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, maxit = 100, residCovWeighted = TRUE,
   useMatrix = useMatrix )
summary( fitsuri4w )


## *************** iterated SUR with 2 restrictions via R and restrict.regMat ****************
fitsurio5 <- systemfit( system, "SUR", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, maxit = 100,
   useMatrix = useMatrix )
print( summary( fitsurio5 ) )
fitsuri5 <- systemfit( system2, "SUR", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, maxit = 100,
   useMatrix = useMatrix )
print( summary( fitsuri5 ) )

## ********* iterated SUR with 2 restrictions via R and restrict.regMat (EViews-like) **********
fitsurio5e <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "noDfCor",
   restrict.matrix = restr3m, restrict.rhs = restr3q, restrict.regMat = tc,
   maxit = 100, useMatrix = useMatrix )
print( summary( fitsurio5e ) )
fitsuri5e <- systemfit( system2, "SUR", data = Kmenta, methodResidCov = "noDfCor",
   restrict.matrix = restr3m, restrict.rhs = restr3q, restrict.regMat = tc,
   maxit = 100, useMatrix = useMatrix )
print( summary( fitsuri5e ) )
nobs( fitsuri5e )

## ********* iterated SUR with 2 restrictions via R and restrict.regMat (methodResidCov="Theil") **********
fitsurio5r2 <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "Theil",
   restrict.matrix = restr3m, restrict.rhs = restr3q, restrict.regMat = tc,
   maxit = 100, x = TRUE, useMatrix = useMatrix )
print( summary( fitsurio5r2 ) )
fitsuri5r2 <- systemfit( system2, "SUR", data = Kmenta, methodResidCov = "Theil",
   restrict.matrix = restr3m, restrict.rhs = restr3q, restrict.regMat = tc,
   maxit = 100, x = TRUE, useMatrix = useMatrix )
print( summary( fitsuri5r2 ) )

## ********* iterated SUR with 2 restrictions via R and restrict.regMat (methodResidCov="max") **********
# fitsuri5e <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "max",
#    restrict.matrix = restr3m, restrict.rhs = restr3q, restrict.regMat = tc,
#    maxit = 100, useMatrix = useMatrix )
# print( summary( fitsuri5e ) )
# print( round( vcov( fitsuri5e ), digits = 6 ) )
# disabled, because the estimation does not converge

## ********* iterated WSUR with 2 restrictions via R and restrict.regMat (methodResidCov="Theil") **********
fitsurio5wr2 <- systemfit( system, "SUR", data = Kmenta, methodResidCov = "Theil",
   restrict.matrix = restr3m, restrict.rhs = restr3q, restrict.regMat = tc,
   maxit = 100, residCovWeighted = TRUE, useMatrix = useMatrix )
summary( fitsurio5wr2 )
fitsuri5wr2 <- systemfit( system2, "SUR", data = Kmenta, methodResidCov = "Theil",
   restrict.matrix = restr3m, restrict.rhs = restr3q, restrict.regMat = tc,
   maxit = 100, residCovWeighted = TRUE, useMatrix = useMatrix )
summary( fitsuri5wr2 )


## *********** estimations with a single regressor ************
fitsurS1 <- systemfit(
   list( price ~ consump - 1, farmPrice ~ consump + trend ), "SUR",
   data = Kmenta, useMatrix = useMatrix )
print( summary( fitsurS1 ) )
nobs( fitsurS1 )
fitsurS2 <- systemfit(
   list( consump ~ price - 1, consump ~ trend - 1 ), "SUR",
   data = Kmenta, useMatrix = useMatrix )
print( summary( fitsurS2 ) )
nobs( fitsurS2 )
fitsurS3 <- systemfit(
   list( consump ~ trend - 1, price ~ trend - 1 ), "SUR",
   data = Kmenta, useMatrix = useMatrix )
nobs( fitsurS3 )
print( summary( fitsurS3 ) )
fitsurS4 <- systemfit(
   list( consump ~ trend - 1, price ~ trend - 1 ), "SUR",
   data = Kmenta, restrict.matrix = matrix( c( 1, -1 ), nrow = 1 ),
   useMatrix = useMatrix )
print( summary( fitsurS4 ) )
nobs( fitsurS4 )
fitsurS5 <- systemfit(
   list( consump ~ 1, price ~ 1 ), "SUR",
   data = Kmenta, useMatrix = useMatrix )
print( summary( fitsurS5 ) )
nobs( fitsurS5 )


## **************** shorter summaries **********************
print( summary( fitsur1e2, useDfSys = TRUE, equations = FALSE ) )

print( summary( fitsur2e, useDfSys = FALSE, residCov = FALSE ) )

print( summary( fitsur3 ), equations = FALSE )

print( summary( fitsur4r3 ), residCov = FALSE, equations = FALSE )

print( summary( fitsur5, residCov = FALSE ), equations = FALSE )

print( summary( fitsur5w, equations = FALSE, residCov = FALSE ),
   equations = TRUE )

print( summary( fitsuri1r3, useDfSys = FALSE ), residCov = FALSE )

print( summary( fitsuri2 ), residCov = FALSE )

print( summary( fitsuri3e, residCov = FALSE, equations = FALSE ) )

print( summary( fitsurio4, residCov = FALSE ), equations = FALSE )
print( summary( fitsuri4, equations = FALSE ), residCov = FALSE )

print( summary( fitsuri4w, useDfSys = FALSE, equations = FALSE ) )

print( summary( fitsurio5r2, equations = FALSE ) )
print( summary( fitsuri5r2 ), residCov = FALSE )


## ****************** residuals **************************
print( residuals( fitsur1e2 ) )
print( residuals( fitsur1e2$eq[[ 2 ]] ) )

print( residuals( fitsur1w ) )
print( residuals( fitsur1w$eq[[ 2 ]] ) )

print( residuals( fitsur2e ) )
print( residuals( fitsur2e$eq[[ 1 ]] ) )

print( residuals( fitsur3 ) )
print( residuals( fitsur3$eq[[ 2 ]] ) )

print( residuals( fitsur4r3 ) )
print( residuals( fitsur4r3$eq[[ 1 ]] ) )

print( residuals( fitsur5 ) )
print( residuals( fitsur5$eq[[ 2 ]] ) )

print( residuals( fitsuri1r3 ) )
print( residuals( fitsuri1r3$eq[[ 1 ]] ) )

print( residuals( fitsuri2 ) )
print( residuals( fitsuri2$eq[[ 2 ]] ) )

print( residuals( fitsuri3e ) )
print( residuals( fitsuri3e$eq[[ 1 ]] ) )

print( residuals( fitsurio4 ) )
print( residuals( fitsurio4$eq[[ 2 ]] ) )
print( residuals( fitsuri4 ) )
print( residuals( fitsuri4$eq[[ 2 ]] ) )

print( residuals( fitsuri4w ) )
print( residuals( fitsuri4w$eq[[ 2 ]] ) )

print( residuals( fitsurio5r2 ) )
print( residuals( fitsurio5r2$eq[[ 1 ]] ) )
print( residuals( fitsuri5r2 ) )
print( residuals( fitsuri5r2$eq[[ 1 ]] ) )


## *************** coefficients *********************
print( round( coef( fitsur1r3 ), digits = 6 ) )
print( round( coef( fitsur1r3$eq[[ 2 ]] ), digits = 6 ) )

print( round( coef( fitsuri2 ), digits = 6 ) )
print( round( coef( fitsuri2$eq[[ 1 ]] ), digits = 6 ) )

print( round( coef( fitsur2we ), digits = 6 ) )
print( round( coef( fitsur2we$eq[[ 1 ]] ), digits = 6 ) )

print( round( coef( fitsur3 ), digits = 6 ) )
print( round( coef( fitsur3, modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( fitsur3$eq[[ 2 ]] ), digits = 6 ) )

print( round( coef( fitsur4r2 ), digits = 6 ) )
print( round( coef( fitsur4r2$eq[[ 1 ]] ), digits = 6 ) )

print( round( coef( fitsuri5e ), digits = 6 ) )
print( round( coef( fitsuri5e, modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( fitsuri5e$eq[[ 2 ]] ), digits = 6 ) )

print( round( coef( fitsur5w ), digits = 6 ) )
print( round( coef( fitsur5w, modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( fitsur5w$eq[[ 1 ]] ), digits = 6 ) )


## *************** coefficients with stats *********************
print( round( coef( summary( fitsur1r3 ) ), digits = 6 ) )
print( round( coef( summary( fitsur1r3$eq[[ 2 ]] ) ), digits = 6 ) )

print( round( coef( summary( fitsuri2, useDfSys = FALSE ) ), digits = 6 ) )
print( round( coef( summary( fitsuri2$eq[[ 1 ]], useDfSys = FALSE ) ),
   digits = 6 ) )

print( round( coef( summary( fitsur3 ) ), digits = 6 ) )
print( round( coef( summary( fitsur3 ), modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( summary( fitsur3$eq[[ 2 ]] ) ), digits = 6 ) )

print( round( coef( summary( fitsuri3we ) ), digits = 6 ) )
print( round( coef( summary( fitsuri3we ), modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( summary( fitsuri3we$eq[[ 1 ]] ) ), digits = 6 ) )

print( round( coef( summary( fitsur4r2 ) ), digits = 6 ) )
print( round( coef( summary( fitsur4r2$eq[[ 1 ]] ) ), digits = 6 ) )

print( round( coef( summary( fitsur4we ) ), digits = 6 ) )
print( round( coef( summary( fitsur4we$eq[[ 2 ]] ) ), digits = 6 ) )

print( round( coef( summary( fitsuri5e, useDfSys = FALSE ) ), digits = 6 ) )
print( round( coef( summary( fitsuri5e, useDfSys = FALSE ),
   modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( summary( fitsuri5e$eq[[ 2 ]], useDfSys = FALSE ) ),
   digits = 6 ) )


## *********** variance covariance matrix of the coefficients *******
print( round( vcov( fitsur1e2 ), digits = 6 ) )
print( round( vcov( fitsur1e2$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitsur1r3 ), digits = 6 ) )
print( round( vcov( fitsur1r3$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitsur2e ), digits = 6 ) )
print( round( vcov( fitsur2e$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitsur3 ), digits = 6 ) )
print( round( vcov( fitsur3, modified.regMat = TRUE ), digits = 6 ) )
print( round( vcov( fitsur3$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitsur3w ), digits = 6 ) )
print( round( vcov( fitsur3w, modified.regMat = TRUE ), digits = 6 ) )
print( round( vcov( fitsur3w$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitsur4r2 ), digits = 6 ) )
print( round( vcov( fitsur4r2$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitsur5e ), digits = 6 ) )
print( round( vcov( fitsur5e, modified.regMat = TRUE ), digits = 6 ) )
print( round( vcov( fitsur5e$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitsuri1r3 ), digits = 6 ) )
print( round( vcov( fitsuri1r3$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitsuri2 ), digits = 6 ) )
print( round( vcov( fitsuri2$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitsuri3e ), digits = 6 ) )
print( round( vcov( fitsuri3e, modified.regMat = TRUE ), digits = 6 ) )
print( round( vcov( fitsuri3e$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitsurio4e ), digits = 6 ) )
print( round( vcov( fitsurio4e$eq[[ 2 ]] ), digits = 6 ) )
print( round( vcov( fitsuri4e ), digits = 6 ) )
print( round( vcov( fitsuri4e$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitsurio5r2 ), digits = 6 ) )
print( round( vcov( fitsurio5r2, modified.regMat = TRUE ), digits = 6 ) )
print( round( vcov( fitsurio5r2$eq[[ 1 ]] ), digits = 6 ) )
print( round( vcov( fitsuri5r2 ), digits = 6 ) )
print( round( vcov( fitsuri5r2, modified.regMat = TRUE ), digits = 6 ) )
print( round( vcov( fitsuri5r2$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitsurio5wr2 ), digits = 6 ) )
print( round( vcov( fitsurio5wr2, modified.regMat = TRUE ), digits = 6 ) )
print( round( vcov( fitsurio5wr2$eq[[ 2 ]] ), digits = 6 ) )


## *********** confidence intervals of coefficients *************
print( confint( fitsur1e2, useDfSys = TRUE ) )
print( confint( fitsur1e2$eq[[ 2 ]], level = 0.9, useDfSys = TRUE ) )

print( confint( fitsur1we2, useDfSys = TRUE ) )
print( confint( fitsur1we2$eq[[ 1 ]], level = 0.9, useDfSys = TRUE ) )

print( confint( fitsur2e, level = 0.9 ) )
print( confint( fitsur2e$eq[[ 1 ]], level = 0.99 ) )

print( confint( fitsur3, level = 0.99 ) )
print( confint( fitsur3$eq[[ 2 ]], level = 0.5 ) )

print( confint( fitsur4r3, level = 0.5 ) )
print( confint( fitsur4r3$eq[[ 1 ]], level = 0.25 ) )

print( confint( fitsur5, level = 0.25 ) )
print( confint( fitsur5$eq[[ 2 ]], level = 0.975 ) )

print( confint( fitsuri1r3, level = 0.975 ) )
print( confint( fitsuri1r3$eq[[ 1 ]], level = 0.999 ) )

print( confint( fitsuri2, level = 0.999 ) )
print( confint( fitsuri2$eq[[ 2 ]], level = 0.1 ) )

print( confint( fitsuri3e, level = 0.1 ) )
print( confint( fitsuri3e$eq[[ 1 ]], level = 0.01 ) )

print( confint( fitsurio4, level = 0.01 ) )
print( confint( fitsurio4$eq[[ 2 ]], level = 0.33 ) )
print( confint( fitsuri4, level = 0.01 ) )
print( confint( fitsuri4$eq[[ 2 ]], level = 0.33 ) )

print( confint( fitsurio4w, level = 0.01 ) )
print( confint( fitsurio4w$eq[[ 1 ]], level = 0.33 ) )

print( confint( fitsurio5r2, level = 0.33 ) )
print( confint( fitsurio5r2$eq[[ 1 ]] ) )
print( confint( fitsuri5r2, level = 0.33 ) )
print( confint( fitsuri5r2$eq[[ 1 ]] ) )


## *********** fitted values *************
print( fitted( fitsur1e2 ) )
print( fitted( fitsur1e2$eq[[ 2 ]] ) )

print( fitted( fitsur2e ) )
print( fitted( fitsur2e$eq[[ 1 ]] ) )

print( fitted( fitsur2we ) )
print( fitted( fitsur2we$eq[[ 2 ]] ) )

print( fitted( fitsur3 ) )
print( fitted( fitsur3$eq[[ 2 ]] ) )

print( fitted( fitsur4r3 ) )
print( fitted( fitsur4r3$eq[[ 1 ]] ) )

print( fitted( fitsur5 ) )
print( fitted( fitsur5$eq[[ 2 ]] ) )

print( fitted( fitsuri1r3 ) )
print( fitted( fitsuri1r3$eq[[ 1 ]] ) )

print( fitted( fitsuri1wr3 ) )
print( fitted( fitsuri1wr3$eq[[ 2 ]] ) )

print( fitted( fitsuri2 ) )
print( fitted( fitsuri2$eq[[ 2 ]] ) )

print( fitted( fitsuri3e ) )
print( fitted( fitsuri3e$eq[[ 1 ]] ) )

print( fitted( fitsurio4 ) )
print( fitted( fitsurio4$eq[[ 2 ]] ) )
print( fitted( fitsuri4 ) )
print( fitted( fitsuri4$eq[[ 2 ]] ) )

print( fitted( fitsurio5r2 ) )
print( fitted( fitsurio5r2$eq[[ 1 ]] ) )
print( fitted( fitsuri5r2 ) )
print( fitted( fitsuri5r2$eq[[ 1 ]] ) )


## *********** predicted values *************
predictData <- Kmenta
predictData$consump <- NULL
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fitsur1e2, se.fit = TRUE, interval = "prediction",
   useDfSys = TRUE ) )
print( predict( fitsur1e2$eq[[ 2 ]], se.fit = TRUE, interval = "prediction",
   useDfSys = TRUE ) )

print( predict( fitsur2e, se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData ) )
print( predict( fitsur2e$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData ) )

print( predict( fitsur3, se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )
print( predict( fitsur3$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )

print( predict( fitsur4r3, se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )
print( predict( fitsur4r3$eq[[ 1 ]], se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )

print( predict( fitsur4we, se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )
print( predict( fitsur4we$eq[[ 2 ]], se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )

print( predict( fitsur5, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )
print( predict( fitsur5$eq[[ 2 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )

print( predict( fitsuri1r3, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )
print( predict( fitsuri1r3$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )

print( predict( fitsuri2, se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )
print( predict( fitsuri2$eq[[ 2 ]], se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )

print( predict( fitsuri2w, se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )
print( predict( fitsuri2w$eq[[ 2 ]], se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )

print( predict( fitsuri3e, interval = "prediction", level = 0.925 ) )
print( predict( fitsuri3e$eq[[ 1 ]], interval = "prediction", level = 0.925 ) )

print( predict( fitsurio4, interval = "confidence", newdata = predictData ) )
print( predict( fitsurio4$eq[[ 2 ]], interval = "confidence",
   newdata = predictData ) )
print( predict( fitsuri4, interval = "confidence", newdata = predictData ) )
print( predict( fitsuri4$eq[[ 2 ]], interval = "confidence",
   newdata = predictData ) )

print( predict( fitsurio5r2 ) )
print( predict( fitsurio5r2$eq[[ 1 ]] ) )
print( predict( fitsuri5r2 ) )
print( predict( fitsuri5r2$eq[[ 1 ]] ) )

# predict just one observation
smallData <- data.frame( price = 130, income = 150, farmPrice = 120,
   trend = 25 )

print( predict( fitsur1e2, newdata = smallData ) )
print( predict( fitsur1e2$eq[[ 1 ]], newdata = smallData ) )

print( predict( fitsur2e, se.fit = TRUE, level = 0.9,
   newdata = smallData ) )
print( predict( fitsur2e$eq[[ 1 ]], se.pred = TRUE, level = 0.99,
   newdata = smallData ) )

print( predict( fitsur3, interval = "prediction", level = 0.975,
   newdata = smallData ) )
print( predict( fitsur3$eq[[ 1 ]], interval = "confidence", level = 0.8,
   newdata = smallData ) )

print( predict( fitsur4r3, se.fit = TRUE, interval = "confidence",
   level = 0.999, newdata = smallData ) )
print( predict( fitsur4r3$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.75, newdata = smallData ) )

print( predict( fitsur5, se.fit = TRUE, interval = "prediction",
   newdata = smallData ) )
print( predict( fitsur5$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   newdata = smallData ) )

print( predict( fitsurio5r2, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = smallData ) )
print( predict( fitsurio5r2$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, newdata = smallData ) )
print( predict( fitsuri5r2, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = smallData ) )
print( predict( fitsuri5r2$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, newdata = smallData ) )

print( predict( fitsuri5wr2, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = smallData ) )
print( predict( fitsuri5wr2$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, newdata = smallData ) )


## ************ correlation of predicted values ***************
print( correlation.systemfit( fitsur1e2, 2, 1 ) )

print( correlation.systemfit( fitsur2e, 1, 2 ) )

print( correlation.systemfit( fitsur3, 2, 1 ) )

print( correlation.systemfit( fitsur3w, 2, 1 ) )

print( correlation.systemfit( fitsur4r3, 1, 2 ) )

print( correlation.systemfit( fitsur5, 2, 1 ) )

print( correlation.systemfit( fitsuri1r3, 1, 2 ) )

print( correlation.systemfit( fitsuri2, 2, 1 ) )

print( correlation.systemfit( fitsuri2w, 1, 2 ) )

print( correlation.systemfit( fitsuri3e, 1, 2 ) )

print( correlation.systemfit( fitsurio4, 2, 1 ) )
print( correlation.systemfit( fitsuri4, 2, 1 ) )

print( correlation.systemfit( fitsurio5r2, 1, 2 ) )
print( correlation.systemfit( fitsuri5r2, 1, 2 ) )


## ************ Log-Likelihood values ***************
print( logLik( fitsur1e2 ) )
print( logLik( fitsur1e2, residCovDiag = TRUE ) )

print( logLik( fitsur2e ) )
print( logLik( fitsur2e, residCovDiag = TRUE ) )

print( logLik( fitsur3 ) )
print( logLik( fitsur3, residCovDiag = TRUE ) )

print( logLik( fitsur4r3 ) )
print( logLik( fitsur4r3, residCovDiag = TRUE ) )

print( logLik( fitsur5 ) )
print( logLik( fitsur5, residCovDiag = TRUE ) )

print( logLik( fitsur5w ) )
print( logLik( fitsur5w, residCovDiag = TRUE ) )

print( logLik( fitsuri1r3 ) )
print( logLik( fitsuri1r3, residCovDiag = TRUE ) )

print( logLik( fitsuri2 ) )
print( logLik( fitsuri2, residCovDiag = TRUE ) )

print( logLik( fitsuri3e ) )
print( logLik( fitsuri3e, residCovDiag = TRUE ) )

print( logLik( fitsurio4 ) )
print( logLik( fitsurio4, residCovDiag = TRUE ) )

print( logLik( fitsuri4 ) )
print( logLik( fitsuri4, residCovDiag = TRUE ) )

print( logLik( fitsuri4w ) )
print( logLik( fitsuri4w, residCovDiag = TRUE ) )

print( logLik( fitsurio5r2 ) )
print( logLik( fitsurio5r2, residCovDiag = TRUE ) )

print( logLik( fitsuri5r2 ) )
print( logLik( fitsuri5r2, residCovDiag = TRUE ) )


## *********** likelihood ratio tests *************
# testing first restriction
# non-iterating, methodResidCov = 1
print( lrtest( fitsur2, fitsur1 ) )
print( lrtest( fitsur3, fitsur1 ) )
# non-iterating, methodResidCov = 0
print( lrtest( fitsur2e, fitsur1e ) )
print( lrtest( fitsur3e, fitsur1e ) )
# iterating, methodResidCov = 1
print( lrtest( fitsuri2, fitsuri1 ) )
print( lrtest( fitsuri3, fitsuri1 ) )
# iterating, methodResidCov = 0
print( lrtest( fitsuri2e, fitsuri1e ) )
print( lrtest( fitsuri3e, fitsuri1e ) )
# non-iterating, methodResidCov = 1, WSUR
print( lrtest( fitsur3w, fitsur1w ) )

# testing second restriction
# non-iterating, methodResidCov = 1
print( lrtest( fitsur4, fitsur2 ) )
print( lrtest( fitsur4, fitsur3 ) )
print( lrtest( fitsur5, fitsur2 ) )
print( lrtest( fitsur5, fitsur3 ) )
# non-iterating, methodResidCov = 0
print( lrtest( fitsur4e, fitsur2e ) )
print( lrtest( fitsur4e, fitsur3e ) )
print( lrtest( fitsur5e, fitsur2e ) )
print( lrtest( fitsur5e, fitsur3e ) )
# iterating, methodResidCov = 1
print( lrtest( fitsurio4, fitsuri2 ) )
print( lrtest( fitsurio4, fitsuri3 ) )
print( lrtest( fitsurio5, fitsuri2 ) )
print( lrtest( fitsurio5, fitsuri3 ) )
   # corrected
print( lrtest( fitsuri2, fitsuri4 ) )
print( lrtest( fitsuri3, fitsuri4 ) )
print( lrtest( fitsuri2, fitsuri5 ) )
print( lrtest( fitsuri3, fitsuri5 ) )

# iterating, methodResidCov = 0
print( lrtest( fitsurio4e, fitsuri2e ) )
print( lrtest( fitsurio4e, fitsuri3e ) )
print( lrtest( fitsurio5e, fitsuri2e ) )
print( lrtest( fitsurio5e, fitsuri3e ) )
   # corrected
print( lrtest( fitsuri2e, fitsuri4e ) )
print( lrtest( fitsuri3e, fitsuri4e ) )
print( lrtest( fitsuri2e, fitsuri5e ) )
print( lrtest( fitsuri3e, fitsuri5e ) )

# non-iterating, methodResidCov = 0, WSUR
print( lrtest( fitsur4we, fitsur2we ) )

# iterating, methodResidCov = 1, WSUR
print( lrtest( fitsuri2w, fitsuri4w ) )

# testing both of the restrictions
# non-iterating, methodResidCov = 1
print( lrtest( fitsur4, fitsur1 ) )
print( lrtest( fitsur5, fitsur1 ) )
# non-iterating, methodResidCov = 0
print( lrtest( fitsur4e, fitsur1e ) )
print( lrtest( fitsur5e, fitsur1e ) )
# iterating, methodResidCov = 1
print( lrtest( fitsurio4, fitsuri1 ) )
print( lrtest( fitsurio5, fitsuri1 ) )
   # corrected
print( lrtest( fitsuri1, fitsuri4 ) )
print( lrtest( fitsuri1, fitsuri5 ) )
# iterating, methodResidCov = 0
print( lrtest( fitsurio4e, fitsuri1e ) )
print( lrtest( fitsurio5e, fitsuri1e ) )
   # corrected
print( lrtest( fitsuri1e, fitsuri4e ) )
print( lrtest( fitsuri1e, fitsuri5e ) )
# non-iterating, methodResidCov = 1, WSUR
print( lrtest( fitsur5w, fitsur1w ) )

# testing the two restrictions with one call
# non-iterating, methodResidCov = 1
print( lrtest( fitsur4, fitsur2, fitsur1 ) )
print( lrtest( fitsur5, fitsur3, fitsur1 ) )
print( lrtest( fitsur1, fitsur3, fitsur5 ) )
print( lrtest( object = fitsur5, fitsur3, fitsur1 ) )
print( lrtest( fitsur3, object = fitsur5, fitsur1 ) )
print( lrtest( fitsur3, fitsur1, object = fitsur5 ) )
# iterating, methodResidCov = 0
print( lrtest( fitsuri4e, fitsuri2e, fitsuri1e ) )
print( lrtest( fitsuri5e, fitsuri3e, fitsuri1e ) )

## ************** F tests ****************
# testing first restriction
print( linearHypothesis( fitsur1, restrm ) )
linearHypothesis( fitsur1, restrict )

print( linearHypothesis( fitsur1r2, restrm ) )
linearHypothesis( fitsur1r2, restrict )

print( linearHypothesis( fitsuri1e2, restrm ) )
linearHypothesis( fitsuri1e2, restrict )

print( linearHypothesis( fitsuri1r3, restrm ) )
linearHypothesis( fitsuri1r3, restrict )

print( linearHypothesis( fitsur1we2, restrm ) )
linearHypothesis( fitsur1we2, restrict )

print( linearHypothesis( fitsuri1wr3, restrm ) )
linearHypothesis( fitsuri1wr3, restrict )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
restrictOnly2 <- "- demand_price + supply_price = 0.5"
restrictOnly2i <- "- demand_price + supply_income = 0.5"
# first restriction not imposed 
print( linearHypothesis( fitsur1e2, restrOnly2m, restrOnly2q ) )
linearHypothesis( fitsur1e2, restrictOnly2 )

print( linearHypothesis( fitsuri1, restrOnly2m, restrOnly2q ) )
linearHypothesis( fitsuri1, restrictOnly2i )

# first restriction imposed
print( linearHypothesis( fitsur2, restrOnly2m, restrOnly2q ) )
linearHypothesis( fitsur2, restrictOnly2 )

print( linearHypothesis( fitsur3, restrOnly2m, restrOnly2q ) )
linearHypothesis( fitsur3, restrictOnly2 )

print( linearHypothesis( fitsuri2e, restrOnly2m, restrOnly2q ) )
linearHypothesis( fitsuri2e, restrictOnly2i )

print( linearHypothesis( fitsuri3e, restrOnly2m, restrOnly2q ) )
linearHypothesis( fitsuri3e, restrictOnly2i )

print( linearHypothesis( fitsur2we, restrOnly2m, restrOnly2q ) )
linearHypothesis( fitsur2we, restrictOnly2 )

print( linearHypothesis( fitsuri3we, restrOnly2m, restrOnly2q ) )
linearHypothesis( fitsuri3we, restrictOnly2i )

# testing both of the restrictions
print( linearHypothesis( fitsur1r3, restr2m, restr2q ) )
linearHypothesis( fitsur1r3, restrict2 )

print( linearHypothesis( fitsuri1e2, restr2m, restr2q ) )
linearHypothesis( fitsuri1e2, restrict2i )

print( linearHypothesis( fitsur1w, restr2m, restr2q ) )
linearHypothesis( fitsur1w, restrict2 )

print( linearHypothesis( fitsuri1wr3, restr2m, restr2q ) )
linearHypothesis( fitsuri1wr3, restrict2i )


## ************** Wald tests ****************
# testing first restriction
print( linearHypothesis( fitsur1, restrm, test = "Chisq" ) )
linearHypothesis( fitsur1, restrict, test = "Chisq" )

print( linearHypothesis( fitsur1r2, restrm, test = "Chisq" ) )
linearHypothesis( fitsur1r2, restrict, test = "Chisq" )

print( linearHypothesis( fitsuri1e2, restrm, test = "Chisq" ) )
linearHypothesis( fitsuri1e2, restrict, test = "Chisq" )

print( linearHypothesis( fitsuri1r3, restrm, test = "Chisq" ) )
linearHypothesis( fitsuri1r3, restrict, test = "Chisq" )

print( linearHypothesis( fitsur1w, restrm, test = "Chisq" ) )
linearHypothesis( fitsur1w, restrict, test = "Chisq" )

# testing second restriction
# first restriction not imposed
print( linearHypothesis( fitsur1e2, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fitsur1e2, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fitsuri1, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fitsuri1, restrictOnly2i, test = "Chisq" )

# first restriction imposed
print( linearHypothesis( fitsur2, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fitsur2, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fitsur3, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fitsur3, restrictOnly2, test = "Chisq" )

print( linearHypothesis( fitsuri2e, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fitsuri2e, restrictOnly2i, test = "Chisq" )

print( linearHypothesis( fitsuri3e, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fitsuri3e, restrictOnly2i, test = "Chisq" )

print( linearHypothesis( fitsuri2w, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fitsuri2w, restrictOnly2i, test = "Chisq" )

print( linearHypothesis( fitsur3w, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linearHypothesis( fitsur3w, restrictOnly2, test = "Chisq" )


# testing both of the restrictions
print( linearHypothesis( fitsur1r3, restr2m, restr2q, test = "Chisq" ) )
linearHypothesis( fitsur1r3, restrict2, test = "Chisq" )

print( linearHypothesis( fitsuri1e2, restr2m, restr2q, test = "Chisq" ) )
linearHypothesis( fitsuri1e2, restrict2i, test = "Chisq" )

print( linearHypothesis( fitsur1we2, restr2m, restr2q, test = "Chisq" ) )
linearHypothesis( fitsur1we2, restrict2, test = "Chisq" )

print( linearHypothesis( fitsuri1wr3, restr2m, restr2q, test = "Chisq" ) )
linearHypothesis( fitsuri1wr3, restrict2i, test = "Chisq" )


## ****************** model frame **************************
print( mf <- model.frame( fitsur1e2 ) )
print( mf1 <- model.frame( fitsur1e2$eq[[ 1 ]] ) )
print( attributes( mf1 )$terms )
print( mf2 <- model.frame( fitsur1e2$eq[[ 2 ]] ) )
print( attributes( mf2 )$terms )

print( all.equal( mf, model.frame( fitsur1w ) ) )
print( all.equal( mf1, model.frame( fitsur1w$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitsur2e ) ) )
print( all.equal( mf1, model.frame( fitsur2e$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitsur3 ) ) )
print( all.equal( mf2, model.frame( fitsur3$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitsur4r3 ) ) )
print( all.equal( mf1, model.frame( fitsur4r3$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitsur4we ) ) )
print( all.equal( mf2, model.frame( fitsur4we$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitsur5 ) ) )
print( all.equal( mf2, model.frame( fitsur5$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitsuri1r3 ) ) )
print( all.equal( mf1, model.frame( fitsuri1r3$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitsuri2 ) ) )
print( all.equal( mf1, model.frame( fitsuri2$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitsuri3e ) ) )
print( all.equal( mf1, model.frame( fitsuri3e$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitsurio4 ) ) )
print( all.equal( mf2, model.frame( fitsurio4$eq[[ 2 ]] ) ) )
print( all.equal( mf, model.frame( fitsuri4 ) ) )
print( all.equal( mf1, model.frame( fitsuri4$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitsurio5r2 ) ) )
print( all.equal( mf1, model.frame( fitsurio5r2$eq[[ 1 ]] ) ) )
print( all.equal( mf, model.frame( fitsuri5r2 ) ) )
print( all.equal( mf1, model.frame( fitsuri5r2$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitsuri5wr2 ) ) )
print( all.equal( mf1, model.frame( fitsuri5wr2$eq[[ 1 ]] ) ) )


## **************** model matrix ************************
# with x (returnModelMatrix) = TRUE
print( !is.null( fitsur1e2$eq[[ 1 ]]$x ) )
print( mm <- model.matrix( fitsur1e2 ) )
print( mm1 <- model.matrix( fitsur1e2$eq[[ 1 ]] ) )
print( mm2 <- model.matrix( fitsur1e2$eq[[ 2 ]] ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fitsur1r2 ) ) )
print( all.equal( mm1, model.matrix( fitsur1r2$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur1r2$eq[[ 2 ]] ) ) )
print( !is.null( fitsur1r2$eq[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fitsur2e$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fitsur2e ) ) )
print( all.equal( mm1, model.matrix( fitsur2e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur2e$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fitsur2 ) ) )
print( all.equal( mm1, model.matrix( fitsur2$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur2$eq[[ 2 ]] ) ) )
print( !is.null( fitsur2$eq[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fitsur2we$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fitsur2we ) ) )
print( all.equal( mm1, model.matrix( fitsur2we$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur2we$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fitsur2 ) ) )
print( all.equal( mm1, model.matrix( fitsur2$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur2$eq[[ 2 ]] ) ) )
print( !is.null( fitsuri2$eq[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fitsur3e$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fitsur3e ) ) )
print( all.equal( mm1, model.matrix( fitsur3e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur3e$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fitsur3 ) ) )
print( all.equal( mm1, model.matrix( fitsur3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur3$eq[[ 2 ]] ) ) )
print( !is.null( fitsur3$eq[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fitsur3w$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fitsur3w ) ) )
print( all.equal( mm1, model.matrix( fitsur3w$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur3w$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fitsur3 ) ) )
print( all.equal( mm1, model.matrix( fitsur3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur3$eq[[ 2 ]] ) ) )
print( !is.null( fitsuri3$eq[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fitsur4r3$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fitsur4r3 ) ) )
print( all.equal( mm1, model.matrix( fitsur4r3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur4r3$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fitsur4we ) ) )
print( all.equal( mm1, model.matrix( fitsur4we$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur4we$eq[[ 2 ]] ) ) )
print( !is.null( fitsur4we$eq[[ 1 ]]$x ) )

# with x (returnModelMatrix) = TRUE
print( !is.null( fitsurio5r2$eq[[ 1 ]]$x ) )
print( !is.null( fitsur5$eq[[ 1 ]]$x ) )
print( all.equal( mm, model.matrix( fitsurio5r2 ) ) )
print( all.equal( mm1, model.matrix( fitsurio5r2$eq[[ 1 ]] ) ) )
print( all.equal( mm, model.matrix( fitsur5 ) ) )
print( all.equal( mm1, model.matrix( fitsur5$eq[[ 1 ]] ) ) )
#print( all.equal( mm2, model.matrix( fitsuri5r2$eq[[ 2 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fitsurio5 ) ) )
print( all.equal( mm1, model.matrix( fitsurio5$eq[[ 1 ]] ) ) )

# with x (returnModelMatrix) = FALSE
print( all.equal( mm, model.matrix( fitsur5w ) ) )
print( all.equal( mm1, model.matrix( fitsur5w$eq[[ 1 ]] ) ) )
#print( all.equal( mm2, model.matrix( fitsuri5r2$eq[[ 1 ]] ) ) )
print( !is.null( fitsurio5$eq[[ 1 ]]$x ) )
print( !is.null( fitsur5w$eq[[ 1 ]]$x ) )


## **************** formulas ************************
formula( fitsur1e2 )
formula( fitsur1e2$eq[[ 2 ]] )

formula( fitsur2e )
formula( fitsur2e$eq[[ 1 ]] )

formula( fitsur2we )
formula( fitsur2we$eq[[ 1 ]] )

formula( fitsur3 )
formula( fitsur3$eq[[ 2 ]] )

formula( fitsur4r3 )
formula( fitsur4r3$eq[[ 1 ]] )

formula( fitsur5 )
formula( fitsur5$eq[[ 2 ]] )

formula( fitsuri1r3 )
formula( fitsuri1r3$eq[[ 1 ]] )

formula( fitsuri2 )
formula( fitsuri2$eq[[ 2 ]] )

formula( fitsuri3e )
formula( fitsuri3e$eq[[ 1 ]] )

formula( fitsurio4 )
formula( fitsurio4$eq[[ 2 ]] )
formula( fitsuri4 )
formula( fitsuri4$eq[[ 2 ]] )

formula( fitsurio5r2 )
formula( fitsurio5r2$eq[[ 1 ]] )
formula( fitsuri5r2 )
formula( fitsuri5r2$eq[[ 1 ]] )

formula( fitsuri5wr2 )
formula( fitsuri5wr2$eq[[ 1 ]] )


## **************** model terms *******************
terms( fitsur1e2 )
terms( fitsur1e2$eq[[ 2 ]] )

terms( fitsur2e )
terms( fitsur2e$eq[[ 1 ]] )

terms( fitsur3 )
terms( fitsur3$eq[[ 2 ]] )

terms( fitsur3w )
terms( fitsur3w$eq[[ 2 ]] )

terms( fitsur4r3 )
terms( fitsur4r3$eq[[ 1 ]] )

terms( fitsur4we )
terms( fitsur4we$eq[[ 1 ]] )

terms( fitsur5 )
terms( fitsur5$eq[[ 2 ]] )

terms( fitsuri1r3 )
terms( fitsuri1r3$eq[[ 1 ]] )

terms( fitsuri2 )
terms( fitsuri2$eq[[ 2 ]] )

terms( fitsuri3e )
terms( fitsuri3e$eq[[ 1 ]] )

terms( fitsurio4 )
terms( fitsurio4$eq[[ 2 ]] )
terms( fitsuri4 )
terms( fitsuri4$eq[[ 2 ]] )

terms( fitsurio5r2 )
terms( fitsurio5r2$eq[[ 1 ]] )
terms( fitsuri5r2 )
terms( fitsuri5r2$eq[[ 1 ]] )


## **************** estfun ************************
library( "sandwich" )

estfun( fitsur1 )
round( colSums( estfun( fitsur1 ) ), digits = 7 )

estfun( fitsur1e2 )
round( colSums( estfun( fitsur1e2 ) ), digits = 7 )

estfun( fitsur1r3 )
round( colSums( estfun( fitsur1r3 ) ), digits = 7 )

estfun( fitsur1w )
round( colSums( estfun( fitsur1w ) ), digits = 7 )

estfun( fitsuri1e )
round( colSums( estfun( fitsuri1e ) ), digits = 7 )

estfun( fitsuri1wr3 )
round( colSums( estfun( fitsuri1wr3 ) ), digits = 7 )

estfun( fitsurS1 )
round( colSums( estfun( fitsurS1 ) ), digits = 7 )

estfun( fitsurS2 )
round( colSums( estfun( fitsurS2 ) ), digits = 7 )

estfun( fitsurS3 )
round( colSums( estfun( fitsurS3 ) ), digits = 7 )

try( estfun( fitsurS4 ) )

estfun( fitsurS5 )
round( colSums( estfun( fitsurS5 ) ), digits = 7 )


## **************** bread ************************
round( bread( fitsur1 ), digits = 7 )

round( bread( fitsur1e2 ), digits = 7 )

round( bread( fitsur1r3 ), digits = 7 )

round( bread( fitsur1w ), digits = 7 )

round( bread( fitsuri1e ), digits = 7 )

round( bread( fitsuri1wr3 ), digits = 7 )

round( bread( fitsurS1 ), digits = 7 )

round( bread( fitsurS2 ), digits = 7 )

round( bread( fitsurS3 ), digits = 7 )

try( bread( fitsurS4 ) )
