# load the micEconCES package
library( "micEconCES" )

# load the data set
data( "GermanIndustry" )

# remove years 1973 - 1975 because of economic disruptions (see Kemfert 1998)
GermanIndustry <- subset( GermanIndustry, year < 1973 | year > 1975, )

# add a time trend (starting with 0)
GermanIndustry$time <- GermanIndustry$year - 1960

# names of inputs
xNames1 <-  c( "K", "E", "A" )
xNames2 <-  c( "K", "A", "E" )
xNames3 <-  c( "E", "A", "K" )


################# econometric estimation with cesEst ##########################

## Nelder-Mead
cesNm1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "NM", control = list( maxit = 5000 ) )
summary( cesNm1 )
cesNm2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "NM", control = list( maxit = 5000 ) )
summary( cesNm2 )
cesNm3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "NM", control = list( maxit = 5000 ) )
summary( cesNm3 )

## Simulated Annealing
cesSann1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "SANN", control = list( maxit = 2e6 ) )
summary( cesSann1 )
cesSann2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "SANN", control = list( maxit = 2e6 ) )
summary( cesSann2 )
cesSann3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "SANN", control = list( maxit = 2e6 ) )
summary( cesSann3 )

## BFGS
cesBfgs1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgs1 )
cesBfgs2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgs2 )
cesBfgs3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgs3 )

## L-BFGS-B
cesBfgsCon1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "L-BFGS-B", control = list( maxit = 5000 )  )
summary( cesBfgsCon1 )
cesBfgsCon2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "L-BFGS-B", control = list( maxit = 5000 )  )
summary( cesBfgsCon2 )
cesBfgsCon3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "L-BFGS-B", control = list( maxit = 5000 )  )
summary( cesBfgsCon3 )

## Levenberg-Marquardt
cesLm1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLm1 )
cesLm2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLm2 )
cesLm3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLm3 )

## Levenberg-Marquardt, multiplicative error term
cesLm1Me <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   multErr = TRUE, control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLm1Me )
cesLm2Me <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   multErr = TRUE, control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLm2Me )
cesLm3Me <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   multErr = TRUE, control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLm3Me )

## Newton-type
cesNewton1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "Newton", iterlim = 500 )
summary( cesNewton1 )
cesNewton2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "Newton", iterlim = 500 )
summary( cesNewton2 )
cesNewton3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "Newton", iterlim = 500 )
summary( cesNewton3 )

## PORT
cesPort1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "PORT", control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPort1 )
cesPort2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "PORT", control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPort2 )
cesPort3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "PORT", control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPort3 )

## PORT, multiplicative error
cesPort1Me <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "PORT", multErr = TRUE,
   control = list( eval.max = 2000, iter.max = 2000 ) )
summary( cesPort1Me )
cesPort2Me <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "PORT", multErr = TRUE,
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPort2Me )
cesPort3Me <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "PORT", multErr = TRUE,
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPort3Me )

## DE
cesDe1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "DE", control = DEoptim.control( trace = FALSE, NP = 500,
   itermax = 1e4 ) )
summary( cesDe1 )
cesDe2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "DE", control = DEoptim.control( trace = FALSE, NP = 500,
   itermax = 1e4 ) )
summary( cesDe2 )
cesDe3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "DE", control = DEoptim.control( trace = FALSE, NP = 500,
   itermax = 1e4 ) )
summary( cesDe3 )

## nls
cesNls1 <- try( cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   vrs = TRUE, method = "nls" ) )
cesNls2 <- try( cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   vrs = TRUE, method = "nls" ) )
cesNls3 <- try( cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   vrs = TRUE, method = "nls" ) )

## NM - Levenberg-Marquardt
cesNmLm1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   start = coef( cesNm1 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesNmLm1 )
cesNmLm2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   start = coef( cesNm2 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesNmLm2 )
cesNmLm3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   start = coef( cesNm3 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesNmLm3 )

## SANN - Levenberg-Marquardt
cesSannLm1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   start = coef( cesSann1 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesSannLm1 )
cesSannLm2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   start = coef( cesSann2 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesSannLm2 )
cesSannLm3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   start = coef( cesSann3 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesSannLm3 )

## DE - Levenberg-Marquardt
cesDeLm1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   start = coef( cesDe1 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesDeLm1 )
cesDeLm2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   start = coef( cesDe2 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesDeLm2 )
cesDeLm3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   start = coef( cesDe3 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesDeLm3 )

## NM - PORT
cesNmPort1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesNm1 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesNmPort1 )
cesNmPort2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesNm2 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesNmPort2 )
cesNmPort3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesNm3 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesNmPort3 )

## SANN - PORT
cesSannPort1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesSann1 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesSannPort1 )
cesSannPort2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesSann2 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesSannPort2 )
cesSannPort3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesSann3 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesSannPort3 )

## DE - PORT
cesDePort1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesDe1 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesDePort1 )
cesDePort2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesDe2 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesDePort2 )
cesDePort3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesDe3 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesDePort3 )


############# estimation with lambda, rho_1, and rho fixed #####################

# removing technological progress using the lambdas of Kemfert (1998)
# (we can do this, because the model has constant returns to scale)
GermanIndustry$K1 <- GermanIndustry$K * exp( 0.0222 * GermanIndustry$time )
GermanIndustry$E1 <- GermanIndustry$E * exp( 0.0222 * GermanIndustry$time )
GermanIndustry$A1 <- GermanIndustry$A * exp( 0.0222 * GermanIndustry$time )
GermanIndustry$K2 <- GermanIndustry$K * exp( 0.0069 * GermanIndustry$time )
GermanIndustry$E2 <- GermanIndustry$E * exp( 0.0069 * GermanIndustry$time )
GermanIndustry$A2 <- GermanIndustry$A * exp( 0.0069 * GermanIndustry$time )
GermanIndustry$K3 <- GermanIndustry$K * exp( 0.00641 * GermanIndustry$time )
GermanIndustry$E3 <- GermanIndustry$E * exp( 0.00641 * GermanIndustry$time )
GermanIndustry$A3 <- GermanIndustry$A * exp( 0.00641 * GermanIndustry$time )

# names of adjusted inputs
xNames1f <-  c( "K1", "E1", "A1" )
xNames2f <-  c( "K2", "A2", "E2" )
xNames3f <-  c( "E3", "A3", "K3" )

## Nelder-Mead, lambda, rho_1, and rho fixed
cesNmFixed1 <- cesEst( "Y", xNames1f, data = GermanIndustry,
   method = "NM", rho1 = 0.5300, rho = 0.1813,
   control = list( maxit = 5000 ) )
summary( cesNmFixed1 )
cesNmFixed2 <- cesEst( "Y", xNames2f, data = GermanIndustry,
   method = "NM", rho1 = 0.2155, rho = 1.1816,
   control = list( maxit = 5000 ) )
summary( cesNmFixed2 )
cesNmFixed3 <- cesEst( "Y", xNames3f, data = GermanIndustry,
   method = "NM", rho1 = 1.3654, rho = 5.8327,
   control = list( maxit = 5000 ) )
summary( cesNmFixed3 )

## BFGS, lambda, rho_1, and rho fixed
cesBfgsFixed1 <- cesEst( "Y", xNames1f, data = GermanIndustry,
   method = "BFGS", rho1 = 0.5300, rho = 0.1813,
   control = list( maxit = 5000 ) )
summary( cesBfgsFixed1 )
cesBfgsFixed2 <- cesEst( "Y", xNames2f, data = GermanIndustry,
   method = "BFGS", rho1 = 0.2155, rho = 1.1816,
   control = list( maxit = 5000 ) )
summary( cesBfgsFixed2 )
cesBfgsFixed3 <- cesEst( "Y", xNames3f, data = GermanIndustry,
   method = "BFGS", rho1 = 1.3654, rho = 5.8327,
   control = list( maxit = 5000 ) )
summary( cesBfgsFixed3 )

## Levenberg-Marquardt, lambda, rho_1, and rho fixed
cesLmFixed1 <- cesEst( "Y", xNames1f, data = GermanIndustry, 
   rho1 = 0.5300, rho = 0.1813,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmFixed1 )
cesLmFixed2 <- cesEst( "Y", xNames2f, data = GermanIndustry, 
   rho1 = 0.2155, rho = 1.1816,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmFixed2 )
cesLmFixed3 <- cesEst( "Y", xNames3f, data = GermanIndustry, 
   rho1 = 1.3654, rho = 5.8327,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmFixed3 )

## Levenberg-Marquardt, lambda, rho_1, and rho fixed, multiplicative error term
cesLmFixed1Me <- cesEst( "Y", xNames1f, data = GermanIndustry, 
   rho1 = 0.5300, rho = 0.1813, multErr = TRUE,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmFixed1Me )
summary( cesLmFixed1Me, rSquaredLog = FALSE )
cesLmFixed2Me <- cesEst( "Y", xNames2f, data = GermanIndustry, 
   rho1 = 0.2155, rho = 1.1816, multErr = TRUE,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmFixed2Me )
summary( cesLmFixed2Me, rSquaredLog = FALSE )
cesLmFixed3Me <- cesEst( "Y", xNames3f, data = GermanIndustry, 
   rho1 = 1.3654, rho = 5.8327, multErr = TRUE,
   control = nls.lm.control( maxiter = 1024, maxfev = 2000 ) )
summary( cesLmFixed3Me )
summary( cesLmFixed3Me, rSquaredLog = FALSE )

## Newton-type, lambda, rho_1, and rho fixed
cesNewtonFixed1 <- cesEst( "Y", xNames1f, data = GermanIndustry,
   method = "Newton", rho1 = 0.5300, rho = 0.1813, iterlim = 500 )
summary( cesNewtonFixed1 )
cesNewtonFixed2 <- cesEst( "Y", xNames2f, data = GermanIndustry,
   method = "Newton", rho1 = 0.2155, rho = 1.1816, iterlim = 500 )
summary( cesNewtonFixed2 )
cesNewtonFixed3 <- cesEst( "Y", xNames3f, data = GermanIndustry,
   method = "Newton", rho1 = 1.3654, rho = 5.8327, iterlim = 500 )
summary( cesNewtonFixed3 )

## PORT, lambda, rho_1, and rho fixed
cesPortFixed1 <- cesEst( "Y", xNames1f, data = GermanIndustry,
   method = "PORT", rho1 = 0.5300, rho = 0.1813,
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortFixed1 )
cesPortFixed2 <- cesEst( "Y", xNames2f, data = GermanIndustry,
   method = "PORT", rho1 = 0.2155, rho = 1.1816,
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortFixed2 )
cesPortFixed3 <- cesEst( "Y", xNames3f, data = GermanIndustry,
   method = "PORT", rho1 = 1.3654, rho = 5.8327,
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortFixed3 )

# compare RSSs of models with lambda, rho_1, and rho fixed
print( matrix( c( cesNmFixed1$rss, cesBfgsFixed1$rss, cesLmFixed1$rss,
   cesNewtonFixed1$rss, cesPortFixed1$rss ), ncol = 1 ), digits = 16 )
cesFixed1 <- cesLmFixed1
print( matrix( c( cesNmFixed2$rss, cesBfgsFixed2$rss, cesLmFixed2$rss,
   cesNewtonFixed2$rss, cesPortFixed2$rss ), ncol = 1 ), digits = 16 )
cesFixed2 <- cesLmFixed2
print( matrix( c( cesNmFixed3$rss, cesBfgsFixed3$rss, cesLmFixed3$rss,
   cesNewtonFixed3$rss, cesPortFixed3$rss ), ncol = 1 ), digits = 16 )
cesFixed3 <- cesLmFixed3

## check if removing the technical progress worked as expected
Y2Calc <- cesCalc( xNames2f, data = GermanIndustry, 
   coef = coef( cesFixed2 ), nested = TRUE )
all.equal( Y2Calc, fitted( cesFixed2 ) )
Y2TcCalc <- cesCalc( sub( "[123]$", "", xNames2f ), tName = "time", 
   data = GermanIndustry, 
   coef = c( coef( cesFixed2 )[1], lambda = 0.0069, coef( cesFixed2 )[-1] ), 
   nested = TRUE )
all.equal( Y2Calc, Y2TcCalc )


########## Grid Search for Rho_1 and Rho ##############
rhoVec <- c( seq( -1, 1, 0.1 ), seq( 1.2, 4, 0.2 ), seq( 4.4, 14, 0.4 ) )

## BFGS, grid search for rho_1 and rho
cesBfgsGridRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, 
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgsGridRho1 )
plot( cesBfgsGridRho1 )
cesBfgsGridRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, 
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgsGridRho2 )
plot( cesBfgsGridRho2 )
cesBfgsGridRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, 
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgsGridRho3 )
plot( cesBfgsGridRho3 )

# BFGS with grid search estimates as starting values
cesBfgsGridStartRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   start = coef( cesBfgsGridRho1 ),
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgsGridStartRho1 )
cesBfgsGridStartRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   start = coef( cesBfgsGridRho2 ), 
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgsGridStartRho2 )
cesBfgsGridStartRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   start = coef( cesBfgsGridRho3 ), 
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgsGridStartRho3 )

## Levenberg-Marquardt, grid search for rho1 and rho
cesLmGridRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmGridRho1 )
plot( cesLmGridRho1 )
cesLmGridRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmGridRho2 )
plot( cesLmGridRho2 )
cesLmGridRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmGridRho3 )
plot( cesLmGridRho3 )

# LM with grid search estimates as starting values
cesLmGridStartRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   start = coef( cesLmGridRho1 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmGridStartRho1 )
cesLmGridStartRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   start = coef( cesLmGridRho2 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmGridStartRho2 )
cesLmGridStartRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   start = coef( cesLmGridRho3 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmGridStartRho3 )

## Newton-type, grid search for rho_1 and rho
cesNewtonGridRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, 
   method = "Newton", iterlim = 500 )
summary( cesNewtonGridRho1 )
plot( cesNewtonGridRho1 )
cesNewtonGridRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, 
   method = "Newton", iterlim = 500, check.analyticals = FALSE )
summary( cesNewtonGridRho2 )
plot( cesNewtonGridRho2 )
cesNewtonGridRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, 
   method = "Newton", iterlim = 500 )
summary( cesNewtonGridRho3 )
plot( cesNewtonGridRho3 )

# Newton-type with grid search estimates as starting values
cesNewtonGridStartRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   start = coef( cesNewtonGridRho1 ),
   method = "Newton", iterlim = 500, check.analyticals = FALSE )
summary( cesNewtonGridStartRho1 )
cesNewtonGridStartRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   start = coef( cesNewtonGridRho2 ),
   method = "Newton", iterlim = 500, check.analyticals = FALSE )
summary( cesNewtonGridStartRho2 )
cesNewtonGridStartRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   start = coef( cesNewtonGridRho3 ),
   method = "Newton", iterlim = 500, check.analyticals = FALSE )
summary( cesNewtonGridStartRho3 )

## PORT, grid search for rho1 and rho
cesPortGridRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, method = "PORT",
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortGridRho1 )
plot( cesPortGridRho1 )
cesPortGridRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, method = "PORT",
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortGridRho2 )
plot( cesPortGridRho2 )
cesPortGridRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, method = "PORT",
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortGridRho3 )
plot( cesPortGridRho3 )

# PORT with grid search estimates as starting values
cesPortGridStartRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   start = coef( cesPortGridRho1 ),, method = "PORT",
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortGridStartRho1 )
cesPortGridStartRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   start = coef( cesPortGridRho2 ), method = "PORT",
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortGridStartRho2 )
cesPortGridStartRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   start = coef( cesPortGridRho3 ), method = "PORT",
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortGridStartRho3 )


########  estimations for different industrial sectors  ##########

# remove years with missing or incomplete data
GermanIndustry <- subset( GermanIndustry, year >= 1970 & year <= 1988, )

# adjust the time trend so that it again starts with 0
GermanIndustry$time <- GermanIndustry$year - 1970

# rhos for grid search
rhoVec <- c( seq( -1, 0.6, 0.2 ), seq( 0.9, 1.5, 0.3 ),
   seq( 2, 10, 1 ), 12, 15, 20, 30, 50, 100 )

# industries (abbreviations)
indAbbr <- c( "C", "S", "N", "I", "V",  "P", "F" )

# names of inputs
xNames <- list()
xNames[[ 1 ]] <-  c( "K", "E", "A" )
xNames[[ 2 ]] <-  c( "K", "A", "E" )
xNames[[ 3 ]] <-  c( "E", "A", "K" )

# list ("folder") for results
indRes <- list()

# names of estimation methods
metNames <- c( "LM", "PORT", "PORT_Grid", "PORT_Start" )

# table for parameter estimates
tabCoef <- array( NA, dim = c( 9, 7, length( metNames ) ),
   dimnames = list( 
      paste( rep( c( "alpha_", "beta_", "m_" ), 3 ), rep( 1:3, each = 3 ),
         " (", rep( c( "rho_1", "rho", "lambda" ), 3 ), ")", sep = "" ),
      indAbbr, metNames ) )

# table for technological change parameters
tabLambda <- array( NA, dim = c( 7, 3, length( metNames ) ),
   dimnames = list( indAbbr, c(1:3), metNames ) )

# table for R-squared values
tabR2 <- tabLambda

# table for RSS values
tabRss <- tabLambda

# table for economic consistency of LM results
tabConsist <- tabLambda[ , , 1, drop = TRUE ]

#econometric estimation with cesEst
for( indNo in 1:length( indAbbr ) ) {

   # name of industry-specific output
   yIndName <- paste( indAbbr[ indNo ], "Y", sep = "_" )

   # sub-list ("subfolder") for all models of this industrie
   indRes[[ indNo ]] <- list()
   
   for( modNo in 1:3 ) {

      cat( "\n=======================================================\n" )
      cat( "Industry No. ", indNo, ", model No. ", modNo, "\n", sep = "" )
      cat( "=======================================================\n\n" )

      # names of industry-specific inputs
      xIndNames <- paste( indAbbr[ indNo ], xNames[[ modNo ]], sep = "_" )
      
      # sub-sub-list for all estimation results of this model/industrie
      indRes[[ indNo ]][[ modNo ]] <- list()
   
      ## Levenberg-Marquardt
      indRes[[ indNo ]][[ modNo ]]$lm <- cesEst( yIndName, xIndNames, 
         tName = "time", data = GermanIndustry,
         control = nls.lm.control( maxiter = 1024, maxfev = 2000 ) )
      print( tmpSum <- summary( indRes[[ indNo ]][[ modNo ]]$lm ) )
      tmpCoef <- coef( indRes[[ indNo ]][[ modNo ]]$lm )
      tabCoef[ ( 3 * modNo - 2 ):( 3 * modNo ), indNo, "LM" ] <- 
         tmpCoef[ c( "rho_1", "rho", "lambda" ) ]
      tabLambda[ indNo, modNo, "LM" ] <- tmpCoef[ "lambda" ]
      tabR2[ indNo, modNo, "LM" ] <- tmpSum$r.squared
      tabRss[ indNo, modNo, "LM" ] <- tmpSum$rss
      tabConsist[ indNo, modNo ] <- tmpCoef[ "gamma" ] >= 0 &
         tmpCoef[ "delta_1" ] >= 0 & tmpCoef[ "delta_1" ] <= 1 &
         tmpCoef[ "delta" ] >= 0 & tmpCoef[ "delta" ] <= 1 &
         tmpCoef[ "rho_1" ] >= -1 & tmpCoef[ "rho" ] >= -1

      ## PORT
      indRes[[ indNo ]][[ modNo ]]$port <- cesEst( yIndName, xIndNames, 
         tName = "time", data = GermanIndustry, method = "PORT",
         control = list( eval.max = 2000, iter.max = 2000 ) )
      print( tmpSum <- summary( indRes[[ indNo ]][[ modNo ]]$port ) )
      tmpCoef <- coef( indRes[[ indNo ]][[ modNo ]]$port )
      tabCoef[ ( 3 * modNo - 2 ):( 3 * modNo ), indNo, "PORT" ] <- 
         tmpCoef[ c( "rho_1", "rho", "lambda" ) ]
      tabLambda[ indNo, modNo, "PORT" ] <- tmpCoef[ "lambda" ]
      tabR2[ indNo, modNo, "PORT" ] <- tmpSum$r.squared
      tabRss[ indNo, modNo, "PORT" ] <- tmpSum$rss

      # PORT, grid search
      indRes[[ indNo ]][[ modNo ]]$portGrid <- cesEst( yIndName, xIndNames, 
         tName = "time", data = GermanIndustry, method = "PORT",
         rho = rhoVec, rho1 = rhoVec,
         control = list( eval.max = 2000, iter.max = 2000 ) )
      print( tmpSum <- summary( indRes[[ indNo ]][[ modNo ]]$portGrid ) )
      tmpCoef <- coef( indRes[[ indNo ]][[ modNo ]]$portGrid )
      tabCoef[ ( 3 * modNo - 2 ):( 3 * modNo ), indNo, "PORT_Grid" ] <- 
         tmpCoef[ c( "rho_1", "rho", "lambda" ) ]
      tabLambda[ indNo, modNo, "PORT_Grid" ] <- tmpCoef[ "lambda" ]
      tabR2[ indNo, modNo, "PORT_Grid" ] <- tmpSum$r.squared
      tabRss[ indNo, modNo, "PORT_Grid" ] <- tmpSum$rss

      # PORT, grid search for starting values
      indRes[[ indNo ]][[ modNo ]]$portStart <- cesEst( yIndName, xIndNames, 
         tName = "time", data = GermanIndustry, method = "PORT",
         start = coef( indRes[[ indNo ]][[ modNo ]]$portGrid ),
         control = list( eval.max = 2000, iter.max = 2000 ) )
      print( tmpSum <- summary( indRes[[ indNo ]][[ modNo ]]$portStart ) )
      tmpCoef <- coef( indRes[[ indNo ]][[ modNo ]]$portStart )
      tabCoef[ ( 3 * modNo - 2 ):( 3 * modNo ), indNo, "PORT_Start" ] <- 
         tmpCoef[ c( "rho_1", "rho", "lambda" ) ]
      tabLambda[ indNo, modNo, "PORT_Start" ] <- tmpCoef[ "lambda" ]
      tabR2[ indNo, modNo, "PORT_Start" ] <- tmpSum$r.squared
      tabRss[ indNo, modNo, "PORT_Start" ] <- tmpSum$rss
   }
}


##########  tables for presenting estimation results ############

result1 <- matrix( NA, nrow = 24, ncol = 9 )
rownames( result1 ) <- c( "Kemfert (1998)", "fixed",
   "Newton", "BFGS", "L-BFGS-B", "PORT", "LM", 
   "NM", "NM - PORT", "NM - LM", 
   "SANN", "SANN - PORT", "SANN - LM", 
   "DE", "DE - PORT", "DE - LM",
   "Newton grid", "Newton grid start", "BFGS grid", "BFGS grid start", 
   "PORT grid", "PORT grid start", "LM grid", "LM grid start" )
colnames( result1 ) <- c( 
   paste( "$\\", names( coef( cesLm1 ) ), "$", sep = "" ), 
   "c", "RSS", "$R^2$" )
result3 <- result2 <- result1

result1[ "Kemfert (1998)", "$\\lambda$" ] <- 0.0222
result2[ "Kemfert (1998)", "$\\lambda$" ] <- 0.0069
result3[ "Kemfert (1998)", "$\\lambda$" ] <- 0.00641

result1[ "Kemfert (1998)", "$\\rho_1$" ] <- 0.5300
result2[ "Kemfert (1998)", "$\\rho_1$" ] <- 0.2155
result3[ "Kemfert (1998)", "$\\rho_1$" ] <- 1.3654

result1[ "Kemfert (1998)", "$\\rho$" ] <- 0.1813
result2[ "Kemfert (1998)", "$\\rho$" ] <- 1.1816
result3[ "Kemfert (1998)", "$\\rho$" ] <- 5.8327

result1[ "Kemfert (1998)", "$R^2$" ] <- 0.9996
result2[ "Kemfert (1998)", "$R^2$" ] <- 0.786
result3[ "Kemfert (1998)", "$R^2$" ] <- 0.9986

result1[ "fixed", ] <- c( coef( cesFixed1 )[1], 0.0222, 
   coef( cesFixed1 )[-1], cesFixed1$convergence,
   cesFixed1$rss, summary( cesFixed1 )$r.squared )
result2[ "fixed", ] <- c( coef( cesFixed2 )[1], 0.0069,
   coef( cesFixed2 )[-1],cesFixed2$convergence,
   cesFixed2$rss, summary( cesFixed2 )$r.squared )
result3[ "fixed", ] <- c( coef( cesFixed3 )[1], 0.00641,
   coef( cesFixed3 )[-1],cesFixed3$convergence,
   cesFixed3$rss, summary( cesFixed3 )$r.squared )

makeRow <- function( model ) {
   if( is.null( model$multErr ) ) {
      model$multErr <- FALSE
   }
   result <- c( coef( model ), 
      ifelse( is.null( model$convergence ), NA, model$convergence ),
      model$rss, summary( model )$r.squared )
   return( result )
}

result1[ "Newton", ] <- makeRow( cesNewton1 )
result2[ "Newton", ] <- makeRow( cesNewton2 )
result3[ "Newton", ] <- makeRow( cesNewton3 )

result1[ "BFGS", ] <- makeRow( cesBfgs1 )
result2[ "BFGS", ] <- makeRow( cesBfgs2 )
result3[ "BFGS", ] <- makeRow( cesBfgs3 )

result1[ "L-BFGS-B", ] <- makeRow( cesBfgsCon1 )
result2[ "L-BFGS-B", ] <- makeRow( cesBfgsCon2 )
result3[ "L-BFGS-B", ] <- makeRow( cesBfgsCon3 )

result1[ "PORT", ] <- makeRow( cesPort1 )
result2[ "PORT", ] <- makeRow( cesPort2 )
result3[ "PORT", ] <- makeRow( cesPort3 )
   
result1[ "LM", ] <- makeRow( cesLm1 )
result2[ "LM", ] <- makeRow( cesLm2 )
result3[ "LM", ] <- makeRow( cesLm3 )

result1[ "NM", ] <- makeRow( cesNm1 )
result2[ "NM", ] <- makeRow( cesNm2 )
result3[ "NM", ] <- makeRow( cesNm3 )

result1[ "NM - LM", ] <- makeRow( cesNmLm1 )
result2[ "NM - LM", ] <- makeRow( cesNmLm2 )
result3[ "NM - LM", ] <- makeRow( cesNmLm3 )

result1[ "NM - PORT", ] <- makeRow( cesNmPort1 )
result2[ "NM - PORT", ] <- makeRow( cesNmPort2 )
result3[ "NM - PORT", ] <- makeRow( cesNmPort3 )

result1[ "SANN", ] <- makeRow( cesSann1 )
result2[ "SANN", ] <- makeRow( cesSann2 )
result3[ "SANN", ] <- makeRow( cesSann3 )

result1[ "SANN - LM", ] <- makeRow( cesSannLm1 )
result2[ "SANN - LM", ] <- makeRow( cesSannLm2 )
result3[ "SANN - LM", ] <- makeRow( cesSannLm3 )

result1[ "SANN - PORT", ] <- makeRow( cesSannPort1 )
result2[ "SANN - PORT", ] <- makeRow( cesSannPort2 )
result3[ "SANN - PORT", ] <- makeRow( cesSannPort3 )

result1[ "DE", ] <- makeRow( cesDe1 )
result2[ "DE", ] <- makeRow( cesDe2 )
result3[ "DE", ] <- makeRow( cesDe3 )

result1[ "DE - LM", ] <- makeRow( cesDeLm1 )
result2[ "DE - LM", ] <- makeRow( cesDeLm2 )
result3[ "DE - LM", ] <- makeRow( cesDeLm3 )

result1[ "DE - PORT", ] <- makeRow( cesDePort1 )
result2[ "DE - PORT", ] <- makeRow( cesDePort2 )
result3[ "DE - PORT", ] <- makeRow( cesDePort3 )

result1[ "Newton grid", ] <- makeRow( cesNewtonGridRho1 )
result2[ "Newton grid", ] <- makeRow( cesNewtonGridRho2 )
result3[ "Newton grid", ] <- makeRow( cesNewtonGridRho3 )

result1[ "Newton grid start", ] <- makeRow( cesNewtonGridStartRho1 )
result2[ "Newton grid start", ] <- makeRow( cesNewtonGridStartRho2 )
result3[ "Newton grid start", ] <- makeRow( cesNewtonGridStartRho3 )

result1[ "BFGS grid", ] <- makeRow( cesBfgsGridRho1 )
result2[ "BFGS grid", ] <- makeRow( cesBfgsGridRho2 )
result3[ "BFGS grid", ] <- makeRow( cesBfgsGridRho3 )

result1[ "BFGS grid start", ] <- makeRow( cesBfgsGridStartRho1 )
result2[ "BFGS grid start", ] <- makeRow( cesBfgsGridStartRho2 )
result3[ "BFGS grid start", ] <- makeRow( cesBfgsGridStartRho3 )

result1[ "PORT grid", ] <- makeRow( cesPortGridRho1 )
result2[ "PORT grid", ] <- makeRow( cesPortGridRho2 )
result3[ "PORT grid", ] <- makeRow( cesPortGridRho3 )

result1[ "PORT grid start", ] <- makeRow( cesPortGridStartRho1 )
result2[ "PORT grid start", ] <- makeRow( cesPortGridStartRho2 )
result3[ "PORT grid start", ] <- makeRow( cesPortGridStartRho3 )

result1[ "LM grid", ] <- makeRow( cesLmGridRho1 )
result2[ "LM grid", ] <- makeRow( cesLmGridRho2 )
result3[ "LM grid", ] <- makeRow( cesLmGridRho3 )

result1[ "LM grid start", ] <- makeRow( cesLmGridStartRho1 )
result2[ "LM grid start", ] <- makeRow( cesLmGridStartRho2 )
result3[ "LM grid start", ] <- makeRow( cesLmGridStartRho3 )

# create LaTeX tables
library( xtable )
colorRows <- function( result ) {
   rownames( result ) <- paste( 
      ifelse( !is.na( result[ , "$\\delta_1$" ] ) & (
         result[ , "$\\delta_1$" ] < 0 | result[ , "$\\delta_1$" ] >1 | 
         result[ , "$\\delta$" ] < 0 | result[ , "$\\delta$" ] > 1 ) | 
         result[ , "$\\rho_1$" ] < -1 | result[ , "$\\rho$" ] < -1, 
         "MarkThisRow ", "" ),
      rownames( result ), sep = "" )
   return( result )
}

printTable <- function( xTab, fileName ) {
   tempFile <- file()
   print( xTab, file = tempFile,
      floating = FALSE, sanitize.text.function = function(x){x} )
   latexLines <- readLines( tempFile )
   close( tempFile )
   for( i in grep( "MarkThisRow ", latexLines, value = FALSE ) ) {
      latexLines[ i ] <- sub( "MarkThisRow", "\\\\color{red}" ,latexLines[ i ] )
      latexLines[ i ] <- gsub( "&", "& \\\\color{red}" ,latexLines[ i ] )
   }
   writeLines( latexLines, fileName )
   invisible( latexLines )
}
   
result1 <- colorRows( result1 )
xTab1 <- xtable( result1, align = "lrrrrrrrrr", 
   digits = c( 0, rep( 4, 6 ), 0, 0, 4 ) )
printTable( xTab1, fileName = "kemfert1Coef.tex" )

result2 <- colorRows( result2 )
xTab2 <- xtable( result2, align = "lrrrrrrrrr", 
   digits = c( 0, rep( 4, 6 ), 0, 0, 4 ) )
printTable( xTab2, fileName = "kemfert2Coef.tex" )

result3 <- colorRows( result3 )
xTab3 <- xtable( result3, align = "lrrrrrrrrr", 
   digits = c( 0, rep( 4, 6 ), 0, 0, 4 ) )
printTable( xTab3, fileName = "kemfert3Coef.tex" )

