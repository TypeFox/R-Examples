lder2.brm <-
function( u, theta, params )
{

## Calculating the probability of response: ##    
  p <- p.brm(theta, params)
  q <- 1 - p

## Calculating the first and second derivatives: ##  
  pder1 <- pder1.brm(theta, params)
  pder2 <- pder2.brm(theta, params)
 
## Calculating two parts of second derivative: ## 
  lder2.1 <- ( -pder1^2 / p^2 ) + ( pder2 / p )
  lder2.2 <- (  pder1^2 / q^2 ) + ( pder2 / q )

  return( u * lder2.1 - (1 - u) * lder2.2 )
    
} # END lder2.brm FUNCTION
