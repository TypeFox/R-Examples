logLik.grm <-
function( u, theta, params,
          type = c("MLE", "BME"),
          ddist = dnorm, ... ) # ddist and ... are prior distribution stuff
{

# u is the response, and x are the parameters.

  type <- type[1]
  
# Then turn params into a matrix and determine stats:
  params <- rbind(params)

## Calculating the loglikelihood without the Bayesian part: ##  
  p      <- p.grm(theta, params) 
  
  logLik <- log( sel.prm(p, u, length(theta), nrow(params), ncol(params)) )


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
  
} # END logLik.grm FUNCTION