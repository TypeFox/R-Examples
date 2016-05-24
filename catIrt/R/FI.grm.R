FI.grm <-
function( params,                            # parameters over which to calculate
          theta,                             # values/estimates of theta
          type = c("expected", "observed"),  # which information to calculate
          resp = NULL )                      # a response vector/matrix            
{
  
# First, make sure that resp is NULL if type is "expected"
  if( type == "expected" )
    resp <- NULL
  
# Then turn params into a matrix and determine stats::
  params <- rbind(params)

  N <- length(theta)
  K <- ncol(params)


#~~~~~~~~~~~~~~~~~#
# Argument Checks #
#~~~~~~~~~~~~~~~~~#

# Then make sure that the arguments are OK:

## 1 ## (Make sure that resp exists if we are calculating observed information)
  if( is.null(resp) & type == "observed" )
    stop( "need response scalar/vector to calculate observed information" )
    
## 2 ## (Make sure that params, theta, resp are ALL numeric)
  if( mode(params) != "numeric" | mode(theta) != "numeric" )
    stop( "params and thet need to be numeric" )
    
  if( !is.null(resp) & mode(resp) != "numeric" )
    stop( "resp needs to be numeric" )
  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Expected Fisher Information #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Expected Fisher Information: sum[P'^2/P]
  if( type == "expected" ){

## Calculating the probability of response: ##  	
  	p <- p.grm(theta, params)

## Calculating the first derivative: ##  	
  	pder1 <- pder1.grm(theta, params)
  	
  	info <- pder1^2 / p
  	
  	if(N == 1){
  	  info <- colSums(info)
    } else{
      info <- do.call(rbind, lapply(split.data.frame(info, rep(1:N, each = K)), colSums))
    } # END ifelse STATEMENT

  } # END if STATEMENT

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Observed Fisher Information #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  if( type == "observed" ){
  	
  	info <- -lder2.grm(resp, theta, params)
    
  } # END if STATEMENT
 

# If theta is a scalar, item information is a vector and test information is a scalar
# If theta is a vector, item information is a matrix and test information is a vector

  if(N == 1){
  
    i.info <- info
    t.info <- sum(info)
    
  } else{
  	
    i.info <- info
    t.info <- rowSums(i.info)
    
  } # END ifelse STATEMENT
              
  sem <- ifelse(test = signif(t.info) > 0, yes = sqrt( 1 / t.info ), no = NA)
               
  return( list(item = drop(i.info), test = t.info, sem = sem, type = type) )
    
} # END FI.grm FUNCTION
