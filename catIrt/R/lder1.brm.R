# l' = sum[ (u - p)*p'/(p*q) ] #

lder1.brm <-
function( u, theta, params,
          type = c("MLE", "WLE")  ) # WLE gives weighted maximum likelihood score fct
{  

# u is the response, theta is ability, and params are the parameters.

  type  <- type[1]
  N     <- length(theta)
  
## Calculating the probability of response: ##
  p  <- p.brm(theta, params)
  q  <- 1 - p
  pq <- p * q 
  
## Calculating the first and second derivatives: ##
  pder1 <- pder1.brm(theta, params)
  pder2 <- pder2.brm(theta, params)
  
## Calculating lder1 for normal/Warm: ##
  if( type == "MLE" ){
  	
  	lder1 <- ( u - p ) * pder1 / pq
    
  } else{
  	
# Calculating Warm correction:
    if(N == 1){
    	  I <- sum( pder1^2 / pq )
    	} else{
      I <- rowSums( pder1^2 / pq )
    } # END ifelse STATEMENT
    
    H <- ( pder1 * pder2 )  / pq
    
    lder1 <- ( u - p ) * pder1 / pq + H / ( 2 * I )
    
  } # END ifelse STATEMENT
  
## Returning Scalar or Vector of logLik's ##
  if(N == 1){
    return( sum(lder1) )
  } else{
    return( rowSums(lder1) )
  } # END ifelse STATEMENT
  
} # END lder1.brm FUNCTION

