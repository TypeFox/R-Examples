## code based on a script written by Geraldine Henningsen. Thanks!

library( "frontier" )
library( "fdrtool" )  # for function rhalfnorm
library( "MCMCpack" )  # for function rdirichlet
options( digits = 5 )

printAll <- function( x ) {
   for( n in names( x ) ) {
      cat( "$", n, "\n", sep = "" )
      if( n %in% c( "olsParam", "gridParam", "mleParam", "olsStdEr", "mleCov" ) ) {
         print( round( x[[ n ]], 2 ) )
      } else if( n %in% c( "fitted", "resid", "olsResid" ) ) {
         print( round( x[[ n ]], 3 ) )
      } else {
         print( x[[ n ]] )
      }
      cat( "\n" )
   }
   cat( "class\n" )
   print( class( x ) )
}

# seed for pseudo random number generator
set.seed( 200 )

# number of observations
nObs <- 500


####################### Generating 3 Input Quantities X #######################
# input quantities
data <- data.frame(
   X1 = rchisq( nObs, 10 ),
   X2 = rchisq( nObs, 10 ),
   X3 = rchisq( nObs, 10 ) )

# logarithms of input quantities
data$lX1 <- log( data$X1 )
data$lX2 <- log( data$X2 )
data$lX3 <- log( data$X3 )


##################### Generating thetas of the 3 Outputs ######################
# proportions of squared output quantities in squared Euclidean distance
# y1^2 + y2^2 + y3^2 = |y|^2
mSquared <- rdirichlet( nObs, c( 1, 1.5, 2 ) )  # yi^2 / |y|^2 

# ratios of output quantities to Euclidean distance
data$m1 <- sqrt( mSquared[ , 1 ] )  # y1 / |y|
data$m2 <- sqrt( mSquared[ , 2 ] )  # y2 / |y|
data$m3 <- sqrt( mSquared[ , 3 ] )  # y3 / |y|

# thetas
data$theta1 <- acos( data$m1 )
data$theta2 <- acos( data$m2 / sin( data$theta1 ) )
# theta3 should always be zero, hence "cos( theta3 )" should always be 1
all.equal( rep( 1, nObs ),
   data$m3 / ( sin( data$theta1 ) * sin( data$theta2 ) ) )


####################### Generating Inefficiency Term u ########################
# parameter of halfnormal distribution for u
# E(x) = 1/theta and Var(x) = (pi-2)/(2*theta^2)
uTheta <- 5 # E(x) = 0.2, Var(x) = 0.02

# inefficiency term u ~ halfnormal
data$u <- rhalfnorm( nObs, theta = uTheta )


########################## Generating Error Term v ############################
# variance of v
vVar <- 0.07

# General Error term v ~ N( 0, sigma )
data$v <- rnorm( nObs, 0, vVar )


############################## Defining Parameters #############################
A_0  <- -2.6

A_1   <-  0.2
A_2   <-  0.4
A_3   <-  0.3
A_T1  <-  0.07
A_T2  <-  0.01

B_1_1  <- -0.003
B_1_2  <-  0.04
B_1_3  <-  0.02
B_1_T1 <- -0.007
B_1_T2 <-  0.004

B_2_2  <- -0.009
B_2_3  <-  0.008
B_2_T1 <- -0.001
B_2_T2 <- -0.002

B_3_3  <- -0.04
B_3_T1 <-  0.011
B_3_T2 <-  0.005

B_T1_T1 <- -0.001
B_T1_T2 <-  0.001

B_T2_T2 <- -0.002


############### checking monotonicity & elast. of scale #######################
# d|y| / dx >= 0
ela <- cbind(
   A_1 + B_1_1 * data$lX1 + B_1_2 * data$lX2 + B_1_3 * data$lX3 +
   B_1_T1 * data$theta1 + B_1_T2 * data$theta2,
   A_2 + B_1_2 * data$lX1 + B_2_2 * data$lX2 + B_2_3 * data$lX3 +
   B_2_T1 * data$theta1 + B_2_T2 * data$theta2,
   A_3 + B_1_3 * data$lX1 + B_2_3 * data$lX2 + B_3_3 * data$lX3 +
   B_3_T1 * data$theta1 + B_3_T2 * data$theta2 )

# checking monotonicity at all data points
colSums( ela > 0 ) # should be nObs, nObs, nObs if monotonicity is fulfilled

# elasticities of scale
hist( rowSums( ela ), plot = FALSE )
# hist( rowSums( ela ) ) 


#################### Calculating logarithm of output distance #################
# dependent variable of the translog stochastic ray frontier production function
data$lY <-  A_0 + A_T1 * data$theta1 + A_T2 * data$theta2 + 
   0.5 * B_T1_T1 * data$theta1^2 + B_T1_T2 * data$theta1 * data$theta2 +
   0.5 * B_T2_T2 * data$theta2^2 +
   A_1 * data$lX1 + A_2 * data$lX2 + A_3 * data$lX3 +
   0.5 * B_1_1 * data$lX1^2 + B_1_2 * data$lX1 * data$lX2 +
   B_1_3 * data$lX1 * data$lX3 + 0.5 * B_2_2 * data$lX2^2 +
   B_2_3 * data$lX2 * data$lX3 + 0.5 * B_3_3 * data$lX3^2 +
   B_1_T1 * data$theta1 * data$lX1 + B_2_T1 * data$theta1 * data$lX2 +
   B_3_T1 * data$theta1 * data$lX3 + B_1_T2 * data$theta2 * data$lX1 + 
   B_2_T2 * data$theta2 * data$lX2 + B_3_T2 * data$theta2 * data$lX3 - 
   data$u + data$v


########### Calculating Distance and Individual Output Quantities #############
# distance (Euclidian distance of all outputs)
data$Y <- exp( data$lY )

# individual output quantities
data$Y1 <- data$m1 * data$Y
data$Y2 <- data$m2 * data$Y
data$Y3 <- data$m3 * data$Y

# some tests
# distance
data$YTest <- sqrt( data$Y1^2 + data$Y2^2 + data$Y3^2 )
all.equal( data$YTest, data$Y )

# theta_1
data$theta1Test <- acos( data$Y1 / data$Y )
all.equal( data$theta1Test, data$theta1 )

# theta_2
data$theta2Test <- acos( data$Y2 / ( data$Y * sin( data$theta1 ) ) )
all.equal( data$theta2Test, data$theta2 )

# theta_3 should always be zero, hence cos( theta_3 ) should always be 1
data$cosTheta3Test <- data$Y3 / 
   ( data$Y * sin( data$theta1 ) * sin( data$theta2 ) )
all.equal( data$cosTheta3Test, rep( 1, nObs ) )


############################# Estimating Model #################################
# estimation
result <- frontierTranslogRay( yNames= c( "Y1", "Y2", "Y3" ), 
   xNames= c( "X1", "X2", "X3" ), data = data )

print( result, digits = 2 )
round( coef( result ), 3 )
print( summary( result ), digits = 1 )
round( efficiencies( result ), 2 )

printAll( result )

all.equal( data$Y, result$distance, tol = 1e-4 )
all.equal( data$theta1, result$theta_1, tol = 1e-4 )
all.equal( data$theta2, result$theta_2, tol = 1e-4 )

# compPlot( exp( -data$u ), efficiencies( result, asInData = TRUE ) )

