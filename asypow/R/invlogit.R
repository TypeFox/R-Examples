invlogit <- function( x ){

# Inverse logit transform of x -- yields probability

infinity <- 78   # Close to largest number that can be exponentiated

tmp <- pmin( infinity, pmax( -infinity, x ) )

tmp <- exp(tmp)

return( tmp/(1+tmp) )     }
