KL.brm <-
function( params,        # parameters over which to calculate
          theta,         # value of theta
          delta = .1 )   # the indifference region specification
{
  
# Turn params into a matrix:
  params <- rbind(params)
  
# Then find the number of people:
  N <- length(theta)
     
#~~~~~~~~~~~~~~~~~#
# Argument Checks #
#~~~~~~~~~~~~~~~~~#
## 1 ## (Make sure that params, thet, resp are ALL numeric)
  if( mode(params) != "numeric" | mode(theta) != "numeric" )
    stop( "params and theta need to be numeric" )
    
#~~~~~~~~~~~~~~~~#
# KL Information #
#~~~~~~~~~~~~~~~~#

## Calculating the prob given particular thetas ##
  p0 <- p.brm(theta - delta, params)
  p1 <- p.brm(theta + delta, params)
  q0 <- 1 - p0
  q1 <- 1 - p1
  
# To prevent computation problems, work with logs and not probabilities:
  info <- p1 * ( log(p1) - log(p0) ) + q1 * ( log(q1) - log(q0) )

  
# If theta is a scalar, item information is a vector and test information is a scalar
# If theta is a vector, item information is a matrix and test information is a vector
  
  if( length(theta) == 1 ){
  
    i.info <- info
    t.info <- sum(info)
    
  } else{
  	
    i.info <- info
    t.info <- rowSums(i.info)
    
  } # END ifelse STATEMENT
                  
                  
  return( list(item = drop(i.info), test = t.info) )
  
} # END KL.brm FUNCTION

