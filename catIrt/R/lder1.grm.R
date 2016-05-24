lder1.grm <-
function( u, theta, params,
          type = c("MLE", "WLE")  ) # WLE gives weighted maximum likelihood score fct
{  

# u is the response, theta is ability, and params are the parameters.

  type  <- type[1]
  
# Then turn params into a matrix and determine stats:
  params <- rbind(params)
  
  N <- length(theta)
  J <- nrow(params)
  K <- ncol(params)
  
## Calculating the probability of response: ##
  p     <- p.grm(theta, params)
  
## Calculating the first and second derivatives: ##
  pder1 <- pder1.grm(theta, params)
  pder2 <- pder2.grm(theta, params)
  
## Calculating lder1 for normal/Warm: ##
  lder1 <- sel.prm(pder1 / p, u, N, J, K)
    
  if( type == "WLE" ){
  	
  	I <- pder1^2 / p
  	H <- pder1 * pder2 / p
  	
# Calculating Warm correction:
    if(N == 1){
      I <- sum(I)
      H <- sum(H)
    } else{
      I <- unlist(lapply(split.data.frame(I, rep(1:N, each = K)), sum))
      H <- unlist(lapply(split.data.frame(H, rep(1:N, each = K)), sum))
    } # END ifelse STATEMENT
    
    lder1 <- lder1 + H / ( 2 * I ) / J
    
  } # END ifelse STATEMENT
  
# Note: The "J" is because lder1 is a vec/mat of length J or ncol J,
#       but H / (2I) is a scal/vec, and we want only ONE H / (2I)
#       to be added to lder1 for each person.

# Note: unlist(lapply) seems to be a bit faster than vapply and a lot
#       faster than laply (from plyr).
  
## Returning Scalar or Vector of logLik's ##
  if( N == 1 ){
    return( sum(lder1) )
  } else{
    return( rowSums(lder1) )
  } # END ifelse STATEMENT
  
} # END lder1.grm FUNCTION

