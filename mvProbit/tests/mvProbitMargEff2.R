# I thank Mohit Batham for providing this R script that demonstrated
# a bug in mvProbitMargEff() when called with 2 dependent variables
library( "mvProbit" )
nObs <- 100
set.seed( 123 )
xData <- data.frame(
  const = rep( 1, nObs ),
  x1 = as.numeric( rnorm( nObs ) > 0 ),
  x2 = as.numeric( rnorm( nObs ) > 0 ),
  x3 = rnorm( nObs ),
  x4 = rnorm( nObs ))
beta <- c( 0.8, 1.2, -1.0, 1.4, -0.8,
  -0.6, 1.0, 0.6, -1.2, -1.6)
sigma <- symMatrix( c( 1, 0.2, 1))
margEffCond <- try( mvProbitMargEff( ~ x1 + x2 + x3 + x4 , coef = beta,
  sigma = sigma, data = xData, cond = TRUE ) )
round( margEffCond, 3 )
