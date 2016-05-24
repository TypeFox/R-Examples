# Examples of using the bootLR package

# Try several different values
testParameters <- matrix( c(
   60, 100, 60, 100,
   99, 100, 50, 100,
   499, 500, 10, 20,
   100, 100, 50, 50,
   20, 20, 350, 500
), nrow=4 )

apply( testParameters, 2, function(x) {
 cat('Running model for:',x,"\n")
  BayesianLR.test( truePos=x[1], totalPos=x[2], trueNeg=x[3], totalNeg=x[4], verbose=TRUE )
} )



# Use returned values for something else

blrt <- BayesianLR.test( 49, 100, 55, 100 )
str(blrt)
cat( "The positive likelihood ratio was", blrt$posLR, "with a confidence interval from", blrt$posLR.ci[1], "to", blrt$posLR.ci[2], ".\n" )


# Examine the stability of the bootstrapped results

lrciStability <- apply( testParameters, 2, function( x, n) {
  cat('Running model for:',x,"\n")
  replicate( n=n, {
    blrt <- BayesianLR.test( truePos=x[1], totalPos=x[2], trueNeg=x[3], totalNeg=x[4], verbose=TRUE )
    with( blrt, c( posLR.ci, negLR.ci ) )
  })
}, n=7 ) # n is number of times to run each set of parameters to examine stability of results