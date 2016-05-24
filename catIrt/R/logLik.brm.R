# l = sum[ u*log(p) + (1 - u)*log(1 - p) ] #

logLik.brm <-
function( u, theta, params,
          type  = c("MLE", "BME"),
          ddist = dnorm, ... ) # ddist and ... are prior distribution stuff
{
  
# u is the response, theta is ability, and params are the parameters.

  type <- type[1]

## Calculating the loglikelihood without the Bayesian part: ##  
  p <- p.brm(theta, params)
  
  if( is.null( dim(u) ) & (length(theta) > 1) & !is.null( dim(params) ) )
    logLik <- t( u * t( log(p) ) + (1 - u) * t( log(1 - p) ) )
  else
    logLik <- u * log(p) + (1 - u) * log(1 - p)


## Now, the Bayesian part: ##  
  if( type == "MLE" )
    bme <- 1
  if( type == "BME" )
    bme <- ddist(x = theta, ... )
  
# if there is a silly prior, set it to something very small
  bme <- ifelse(test = bme <= 0, yes = bme <- 1e-15 , no = bme)
  
## Returning Scalar or Vector of logLik's ##
  if( length(theta) == 1 ){
    return( sum(logLik) + log(bme) )
  } else{
    return( rowSums(logLik) + log(bme) )
  } # END ifelse STATEMENT
  
} # END logLik.brm FUNCTION