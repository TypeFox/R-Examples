### R code from vignette source 'CES.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: CES.Rnw:126-127
###################################################
options( prompt = "R> ", ctinue = "+  ", width = 70, useFancyQuotes = FALSE )


###################################################
### code chunk number 2: CES.Rnw:463-464
###################################################
library( "micEconCES" )


###################################################
### code chunk number 3: CES.Rnw:472-486
###################################################
set.seed( 123 )
cesData <- data.frame(x1 = rchisq(200, 10), x2 = rchisq(200, 10),
   x3 = rchisq(200, 10), x4 = rchisq(200, 10) )
cesData$y2 <- cesCalc( xNames = c( "x1", "x2" ), data = cesData,
   coef = c( gamma = 1, delta = 0.6, rho = 0.5, nu = 1.1 ) )
cesData$y2 <- cesData$y2 + 2.5 * rnorm( 200 )
cesData$y3 <- cesCalc(xNames = c("x1", "x2", "x3"), data = cesData,
   coef = c( gamma = 1, delta_1 = 0.7, delta = 0.6, rho_1 = 0.3, rho = 0.5,
      nu = 1.1), nested = TRUE )
cesData$y3 <- cesData$y3 + 1.5 * rnorm(200)
cesData$y4 <- cesCalc(xNames = c("x1", "x2", "x3", "x4"), data = cesData,
   coef = c(gamma = 1, delta_1 = 0.7, delta_2 = 0.6, delta = 0.5,
   rho_1 = 0.3, rho_2 = 0.4, rho = 0.5, nu = 1.1), nested = TRUE )
cesData$y4 <- cesData$y4 + 1.5 * rnorm(200)


###################################################
### code chunk number 4: CES.Rnw:518-521
###################################################
cesNls <- nls( y2 ~ gamma * ( delta * x1^(-rho) + (1 - delta) * x2^(-rho) )^(-phi / rho),
   data = cesData, start = c( gamma = 0.5, delta = 0.5, rho = 0.25, phi = 1 ) )
print( cesNls )


###################################################
### code chunk number 5: CES.Rnw:646-648
###################################################
cesKmenta <- cesEst( yName = "y2", xNames = c( "x1", "x2" ), data = cesData, 
   method = "Kmenta", vrs = TRUE )


###################################################
### code chunk number 6: CES.Rnw:652-653
###################################################
summary( cesKmenta )


###################################################
### code chunk number 7: CES.Rnw:667-668
###################################################
coef( summary( cesKmenta$kmenta ) )


###################################################
### code chunk number 8: CES.Rnw:680-682
###################################################
pdf( "plotFittedKmenta.pdf", width = 4, height = 4 )
par( mar = c( 4.5, 4, 1, 1 ) )


###################################################
### code chunk number 9: CES.Rnw:684-687
###################################################
library( "miscTools" )
compPlot ( cesData$y2, fitted( cesKmenta ), xlab = "actual values",
   ylab = "fitted values" )


###################################################
### code chunk number 10: CES.Rnw:689-690
###################################################
dev.off()


###################################################
### code chunk number 11: CES.Rnw:766-768
###################################################
cesLm2 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE )
summary( cesLm2 )


###################################################
### code chunk number 12: CES.Rnw:771-773
###################################################
cesLm3 <- cesEst( "y3", c( "x1", "x2", "x3" ), cesData, vrs = TRUE )
summary( cesLm3 )


###################################################
### code chunk number 13: CES.Rnw:776-778
###################################################
cesLm4 <- cesEst( "y4", c( "x1", "x2", "x3", "x4" ), cesData, vrs = TRUE )
summary( cesLm4 )


###################################################
### code chunk number 14: CES.Rnw:784-786
###################################################
pdf( "plotFittedLm.pdf", width = 9, height = 3 )
par( mar = c( 4.5, 4, 3, 1 ), mfrow = c( 1, 3 ) )


###################################################
### code chunk number 15: CES.Rnw:788-794
###################################################
compPlot ( cesData$y2, fitted( cesLm2 ), xlab = "actual values",
   ylab = "fitted values", main = "two-input CES" )
compPlot ( cesData$y3, fitted( cesLm3 ), xlab = "actual values",
   ylab = "fitted values", main = "three-input nested CES" )
compPlot ( cesData$y4, fitted( cesLm4 ), xlab = "actual values",
   ylab = "fitted values", main = "four-input nested CES" )


###################################################
### code chunk number 16: CES.Rnw:796-797
###################################################
dev.off()


###################################################
### code chunk number 17: CES.Rnw:853-855
###################################################
cesCg <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "CG" )
summary( cesCg )


###################################################
### code chunk number 18: CES.Rnw:864-867
###################################################
cesCg2 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "CG",
   control = list( maxit = 1000, reltol = 1e-5 ) )
summary( cesCg2 )


###################################################
### code chunk number 19: CES.Rnw:892-895
###################################################
cesNewton <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE,
   method = "Newton" )
summary( cesNewton )


###################################################
### code chunk number 20: CES.Rnw:925-927
###################################################
cesBfgs <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "BFGS" )
summary( cesBfgs )


###################################################
### code chunk number 21: CES.Rnw:973-976
###################################################
cesNm <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE,
   method = "NM" )
summary( cesNm )


###################################################
### code chunk number 22: CES.Rnw:1004-1006
###################################################
cesSann <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "SANN" )
summary( cesSann )


###################################################
### code chunk number 23: CES.Rnw:1017-1019
###################################################
cesSann2 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "SANN" )
all.equal( cesSann, cesSann2 )


###################################################
### code chunk number 24: CES.Rnw:1024-1033
###################################################
cesSann3 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "SANN",
   random.seed = 1234 )
cesSann4 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "SANN",
   random.seed = 12345 )
cesSann5 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "SANN",
   random.seed = 123456 )
m <- rbind( cesSann = coef( cesSann ), cesSann3 = coef( cesSann3 ),
   cesSann4 = coef( cesSann4 ), cesSann5 = coef( cesSann5 ) )
rbind( m, stdDev = sd( m ) )


###################################################
### code chunk number 25: CES.Rnw:1039-1050
###################################################
cesSannB <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "SANN",
   control = list( maxit = 100000 ) )
cesSannB3 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "SANN",
   random.seed = 1234, control = list( maxit = 100000 ) )
cesSannB4 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "SANN",
   random.seed = 12345, control = list( maxit = 100000 ) )
cesSannB5 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "SANN",
   random.seed = 123456, control = list( maxit = 100000 ) )
m <- rbind( cesSannB = coef( cesSannB ), cesSannB3 = coef( cesSannB3 ),
   cesSannB4 = coef( cesSannB4 ), cesSannB5 = coef( cesSannB5 ) )
rbind( m, stdDev = sd( m ) )


###################################################
### code chunk number 26: CES.Rnw:1080-1083
###################################################
cesDe <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "DE",
   control = list( trace = FALSE ) )
summary( cesDe )


###################################################
### code chunk number 27: CES.Rnw:1090-1093
###################################################
cesDe2 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "DE",
   control = list( trace = FALSE ) )
all.equal( cesDe, cesDe2 )


###################################################
### code chunk number 28: CES.Rnw:1098-1107
###################################################
cesDe3 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "DE",
   random.seed = 1234, control = list( trace = FALSE ) )
cesDe4 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "DE",
   random.seed = 12345, control = list( trace = FALSE ) )
cesDe5 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "DE",
   random.seed = 123456, control = list( trace = FALSE ) )
m <- rbind( cesDe = coef( cesDe ), cesDe3 = coef( cesDe3 ),
   cesDe4 = coef( cesDe4 ), cesDe5 = coef( cesDe5 ) )
rbind( m, stdDev = sd( m ) )


###################################################
### code chunk number 29: CES.Rnw:1119-1129
###################################################
cesDeB <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "DE",
   control = list( trace = FALSE, itermax = 1000 ) )
cesDeB3 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "DE",
   random.seed = 1234, control = list( trace = FALSE, itermax = 1000 ) )
cesDeB4 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "DE",
   random.seed = 12345, control = list( trace = FALSE, itermax = 1000 ) )
cesDeB5 <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE, method = "DE",
   random.seed = 123456, control = list( trace = FALSE, itermax = 1000 ) )
rbind( cesDeB = coef( cesDeB ), cesDeB3 = coef( cesDeB3 ),
   cesDeB4 = coef( cesDeB4 ), cesDeB5 = coef( cesDeB5 ) )


###################################################
### code chunk number 30: CES.Rnw:1182-1185
###################################################
cesLbfgsb <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE,
   method = "L-BFGS-B" )
summary( cesLbfgsb )


###################################################
### code chunk number 31: CES.Rnw:1204-1207
###################################################
cesPort <- cesEst( "y2", c( "x1", "x2" ), cesData, vrs = TRUE,
   method = "PORT" )
summary( cesPort )


###################################################
### code chunk number 32: CES.Rnw:1261-1268
###################################################
cesData$t <- c( 1:200 )
cesData$yt <- cesCalc( xNames = c( "x1", "x2" ), data = cesData, tName = "t",
   coef = c( gamma = 1, delta = 0.6, rho = 0.5, nu = 1.1, lambda = 0.01 ) )
cesData$yt <- cesData$yt + 2.5 * rnorm( 200 ) 
cesTech <- cesEst( "yt", c( "x1", "x2" ), data = cesData, tName = "t",
   vrs = TRUE, method = "LM" )
summary( cesTech )


###################################################
### code chunk number 33: CES.Rnw:1326-1329
###################################################
cesGrid <- cesEst( "y2", c( "x1", "x2" ), data = cesData, vrs = TRUE,
   rho = seq( from = -0.3, to = 1.5, by = 0.1 ) )
summary( cesGrid )


###################################################
### code chunk number 34: CES.Rnw:1339-1341
###################################################
pdf( "plotGrid.pdf", width = 4, height = 4 )
par( mar = c( 4.5, 4, 1, 1 ) )


###################################################
### code chunk number 35: CES.Rnw:1343-1344
###################################################
plot( cesGrid )


###################################################
### code chunk number 36: CES.Rnw:1346-1347
###################################################
dev.off()


###################################################
### code chunk number 37: CES.Rnw:1367-1373
###################################################
ces4Grid <- cesEst( yName = "y4", xNames = c( "x1", "x2", "x3", "x4" ),
   data = cesData, method = "LM",
   rho1 = seq( from = -0.6, to = 0.9, by = 0.3 ),
   rho2 = seq( from = -0.4, to = 0.8, by = 0.2 ),
   rho = seq( from = -0.3, to = 1.7, by = 0.2 ) )
summary( ces4Grid )


###################################################
### code chunk number 38: CES.Rnw:1386-1387
###################################################
pdf("plotGrid4.pdf", width = 4, height = 8 )


###################################################
### code chunk number 39: CES.Rnw:1389-1390
###################################################
plot( ces4Grid )


###################################################
### code chunk number 40: CES.Rnw:1392-1393
###################################################
dev.off()


###################################################
### code chunk number 41: CES.Rnw:1409-1416
###################################################
cesStartGrid <- cesEst( "y2", c( "x1", "x2" ), data = cesData, vrs = TRUE,
   start = coef( cesGrid ) )
summary( cesStartGrid )

ces4StartGrid <- cesEst( "y4", c( "x1", "x2", "x3", "x4" ), data = cesData,
   start = coef( ces4Grid ) )
summary( ces4StartGrid )


###################################################
### code chunk number 42: CES.Rnw:1559-1561
###################################################
pdf( "plotCesCalcRho.pdf", width = 4.2, height = 4 )
par( mar = c( 4.5, 4, 1.5, 1 ) )


###################################################
### code chunk number 43: CES.Rnw:1563-1580
###################################################
rhoData <- data.frame( rho = seq( -2e-6, 2e-6, 5e-9 ),
   yCES = NA, yLin = NA )
# calculate dependent variables
for( i in 1:nrow( rhoData ) ) {
   # vector of coefficients
   cesCoef <- c( gamma = 1, delta = 0.6, rho = rhoData$rho[ i ], nu = 1.1 )
   rhoData$yLin[ i ] <- cesCalc( xNames = c( "x1", "x2" ), data = cesData[1,],
      coef = cesCoef, rhoApprox = Inf )
   rhoData$yCES[ i ] <- cesCalc( xNames = c( "x1", "x2" ), data = cesData[1,],
      coef = cesCoef, rhoApprox = 0 )
}
# normalise output variables
rhoData$yCES <- rhoData$yCES - rhoData$yLin[ rhoData$rho == 0 ]
rhoData$yLin <- rhoData$yLin - rhoData$yLin[ rhoData$rho == 0 ]
plot( rhoData$rho, rhoData$yCES, type = "l", col = "red",
   xlab = "rho", ylab = "y (normalised, red = CES, black = linearised)" )
lines( rhoData$rho, rhoData$yLin )


###################################################
### code chunk number 44: CES.Rnw:1582-1583
###################################################
dev.off()


###################################################
### code chunk number 45: CES.Rnw:1602-1622
###################################################
pdf( "plotCesCalcRho2.pdf", width = 4.2, height = 4 )
par( mar = c( 4.5, 4, 1.5, 1 ) )
rhoData <- data.frame( rho = seq( -1, 3, 5e-2 ),
   yCES = NA, yLin = NA )
# calculate dependent variables
for( i in 1:nrow( rhoData ) ) {
   # vector of coefficients
   cesCoef <- c( gamma = 1, delta = 0.6, rho = rhoData$rho[ i ], nu = 1.1 )
   rhoData$yLin[ i ] <- cesCalc( xNames = c( "x1", "x2" ), data = cesData[1,],
      coef = cesCoef, rhoApprox = Inf )
   rhoData$yCES[ i ] <- cesCalc( xNames = c( "x1", "x2" ), data = cesData[1,],
      coef = cesCoef, rhoApprox = 0 )
}
# normalise output variables
rhoData$yCES <- rhoData$yCES - rhoData$yLin[ rhoData$rho == 0 ]
rhoData$yLin <- rhoData$yLin - rhoData$yLin[ rhoData$rho == 0 ]
plot( rhoData$rho, rhoData$yCES, type = "l", col = "red",
   xlab = "rho", ylab = "y (normalised, red = CES, black = linearised)" )
lines( rhoData$rho, rhoData$yLin )
dev.off()


###################################################
### code chunk number 46: CES.Rnw:2113-2115
###################################################
data( "GrowthDJ", package = "AER" )
GrowthDJ <- subset( GrowthDJ, oil == "no" )


###################################################
### code chunk number 47: CES.Rnw:2122-2124
###################################################
GrowthDJ$x1 <- 1
GrowthDJ$x2 <- ( GrowthDJ$popgrowth + 5 ) / GrowthDJ$invest


###################################################
### code chunk number 48: CES.Rnw:2132-2134
###################################################
cesNls <- cesEst( "gdp85", c( "x1", "x2"), data = GrowthDJ )
summary( cesNls, ela = FALSE)


###################################################
### code chunk number 49: CES.Rnw:2138-2140
###################################################
cat( "alpha =", ( coef( cesNls )[ "delta" ] - 1 ) / coef( cesNls )[ "delta" ], "\n" )
cat( "sigma =", 1 / ( 1 - coef( cesNls )[ "rho" ] ), "\n" )


###################################################
### code chunk number 50: CES.Rnw:2153-2156
###################################################
cdNls <- cesEst( "gdp85", c( "x1", "x2"), data = GrowthDJ, rho = 0 )
summary( cdNls, ela = FALSE )
cat( "alpha =", ( coef( cdNls )[ "delta" ] - 1 ) /  coef( cdNls )[ "delta" ], "\n" )


###################################################
### code chunk number 51: CES.Rnw:2167-2170
###################################################
cdLog <- cesEst( "gdp85", c( "x1", "x2"), data = GrowthDJ, rho = 0, multErr = TRUE )
summary( cdLog, ela = FALSE )
cat( "alpha =", ( coef( cdLog )[ "delta" ] - 1 ) /  coef( cdLog )[ "delta" ], "\n" )


###################################################
### code chunk number 52: CES.Rnw:2249-2252
###################################################
data( "GermanIndustry" )
GermanIndustry$time <- GermanIndustry$year - 1960
GermanIndustry <- subset( GermanIndustry, year < 1973 | year > 1975, )


###################################################
### code chunk number 53: CES.Rnw:2258-2266
###################################################
cesNls1 <- try( nls( Y ~ gamma * exp( lambda * time ) *
   ( delta2 * ( delta1 * K ^(-rho1) +
      ( 1 - delta1 ) * E^(-rho1) )^( rho / rho1 ) +
   ( 1 - delta2 ) * A^(-rho) )^( - 1 / rho ),
   start = c( gamma = 1, lambda = 0.015, delta1 = 0.5,
      delta2 = 0.5, rho1= 0.2, rho = 0.2 ),
   data = GermanIndustry ) )
cat( cesNls1 )


###################################################
### code chunk number 54: CES.Rnw:2274-2278
###################################################
cesLm1 <- cesEst( "Y", c( "K", "E", "A" ), tName = "time",
   data = GermanIndustry,
   control = nls.lm.control( maxiter = 1024, maxfev = 2000 ) )
summary( cesLm1 )


###################################################
### code chunk number 55: CES.Rnw:2318-2322
###################################################
cesPort1 <- cesEst( "Y", c( "K", "E", "A" ), tName = "time",
   data = GermanIndustry, method = "PORT",
   control = list( iter.max = 1000, eval.max = 2000) )
summary( cesPort1 )


###################################################
### code chunk number 56: CES.Rnw:2361-2368
###################################################
GermanIndustry$K1 <- GermanIndustry$K * exp( 0.0222 * GermanIndustry$time )
GermanIndustry$E1 <- GermanIndustry$E * exp( 0.0222 * GermanIndustry$time )
GermanIndustry$A1 <- GermanIndustry$A * exp( 0.0222 * GermanIndustry$time )
cesLmFixed1 <- cesEst( "Y",  c( "K1", "E1", "A1" ), data = GermanIndustry,
   rho1 = 0.5300, rho = 0.1813,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmFixed1 )


###################################################
### code chunk number 57: CES.Rnw:2555-2556
###################################################
rhoVec <- c(seq(-1, 1, 0.1), seq(1.2, 4, 0.2), seq(4.4, 14, 0.4))


###################################################
### code chunk number 58: CES.Rnw:2558-2562 (eval = FALSE)
###################################################
## cesLmGridRho1 <- cesEst("Y", c("K", "E", "A"), tName = "time", 
##    data = GermanIndustry, 
##    control = nls.lm.control(maxiter = 1000, maxfev = 2000), 
##    rho1 = rhoVec, rho = rhoVec)


###################################################
### code chunk number 59: CES.Rnw:2564-2565
###################################################
load(system.file("Kemfert98Nest1GridLm.RData", package = "micEconCES"))


###################################################
### code chunk number 60: CES.Rnw:2567-2568
###################################################
print(cesLmGridRho1)


###################################################
### code chunk number 61: CES.Rnw:2573-2575
###################################################
pdf( "kemfertGrid.pdf", width = 7, height = 7 )
par( mar = c( 0, 1.5, 3, 0 ) )


###################################################
### code chunk number 62: CES.Rnw:2577-2578
###################################################
plot( cesLmGridRho1 )


###################################################
### code chunk number 63: CES.Rnw:2580-2581
###################################################
dev.off()


###################################################
### code chunk number 64: CES.Rnw:2612-2614
###################################################
pdf( "kemfertGridBest.pdf", width = 7, height = 7 )
par( mar = c( 0, 1.5, 3, 0 ) )


###################################################
### code chunk number 65: CES.Rnw:2616-2619
###################################################
cesLmGridRho1a <- cesLmGridRho1
cesLmGridRho1a$rssArray[ cesLmGridRho1a$rssArray >= 5000 ] <- NA
plot( cesLmGridRho1a )


###################################################
### code chunk number 66: CES.Rnw:2621-2622
###################################################
dev.off()


