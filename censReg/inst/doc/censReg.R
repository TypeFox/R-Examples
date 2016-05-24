### R code from vignette source 'censReg.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: censReg.Rnw:51-52
###################################################
options( prompt = "R> ", ctinue = "+  " )


###################################################
### code chunk number 2: censReg.Rnw:213-214
###################################################
library( "censReg" )


###################################################
### code chunk number 3: censReg.Rnw:232-233
###################################################
data( "Affairs", package = "AER" )


###################################################
### code chunk number 4: censReg.Rnw:243-245
###################################################
estResult <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs )


###################################################
### code chunk number 5: censReg.Rnw:250-251
###################################################
summary( estResult )


###################################################
### code chunk number 6: censReg.Rnw:266-268
###################################################
estResultMinus <- censReg( I( - affairs ) ~ age + yearsmarried + religiousness +
   occupation + rating, left = -Inf, right = 0, data = Affairs )


###################################################
### code chunk number 7: censReg.Rnw:273-274
###################################################
cbind( coef( estResult ), coef( estResultMinus ) )


###################################################
### code chunk number 8: censReg.Rnw:498-509
###################################################
set.seed( 123 )
pData <- data.frame(
   id = rep( paste( "F", 1:15, sep = "_" ), each = 4 ),
   time = rep( 1981:1984, 15 ) )
pData$mu <- rep( rnorm( 15 ), each = 4 )
pData$x1 <- rnorm( 60 )
pData$x2 <- runif( 60 )
pData$ys <- -1 + pData$mu + 2 * pData$x1 + 3 * pData$x2 + rnorm( 60 )
pData$y <- ifelse( pData$ys > 0, pData$ys, 0 )
library( plm )
pData <- pdata.frame( pData, c( "id", "time" ) )


###################################################
### code chunk number 9: censReg.Rnw:515-517
###################################################
system.time( panelResult <- censReg( y ~ x1 + x2, data = pData, method = "BHHH" ) )
summary( panelResult )


###################################################
### code chunk number 10: censReg.Rnw:529-538
###################################################
nGHQ <- 2^(2:6)
times <- numeric( length( nGHQ ) )
results <- list()
for( i in 1:length (nGHQ ) ) {
   times[i] <- system.time( results[[i]] <- censReg( y ~ x1 + x2, data = pData,
   method = "BHHH", nGHQ = nGHQ[i] ) )[1]
}
names(results)<-nGHQ
round( rbind(sapply( results, coef ),times),4)


