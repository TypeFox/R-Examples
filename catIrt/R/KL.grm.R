KL.grm <-
function( params,       # parameters over which to calculate
          theta,        # value(s) of theta
          delta = .1)   # the indifference region specification
{
  
# Turn params into a matrix:
  params <- rbind(params)
  
# Then find the number of people:
  N <- length(theta)
  K <- ncol(params)
    
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
  p0 <- p.grm(theta - delta, params)
  p1 <- p.grm(theta + delta, params)
  
# To prevent computation problems, work with logs and not probabilities:
  info <- p1 * ( log(p1) - log(p0) )
 
  if(N == 1){
  	info <- colSums(info)
  } else{
    info <- do.call(rbind, lapply(split.data.frame(info, rep(1:N, each = K)), colSums))
  } # END ifelse STATEMENT
  
# If theta is a scalar, item information is a vector and test information is a scalar
# If theta is a vector, item information is a matrix and test information is a vector
  
  if(N == 1){
  
    i.info <- info
    t.info <- sum(info)
    
  } else{
  	
    i.info <- info
    t.info <- rowSums(i.info)
    
  } # END ifelse STATEMENT
                                 
  return( list(item = drop(i.info), test = t.info) )
  
} # END KL.grm FUNCTION

