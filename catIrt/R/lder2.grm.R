lder2.grm <-
function( u, theta, params )
{

# Then turn params into a matrix and determine stats:
  params <- rbind(params)
  
  N <- length(theta)
  J <- nrow(params)
  K <- ncol(params)

## Calculating the probability of response: ##    
  p <- p.grm(theta, params)

## Calculating the first and second derivatives: ##  
  pder1 <- pder1.grm(theta, params)
  pder2 <- pder2.grm(theta, params)
 
## Calculating two parts of second derivative: ## 
  lder2 <- ( -pder1^2 / p^2 ) + ( pder2 / p )

  return( sel.prm(lder2, u, N, J, K) )
    
} # END lder2.grm FUNCTION

