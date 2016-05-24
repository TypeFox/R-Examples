logit <- function( p ){ 

# Logit transform of probability p

infinity <- 1.0e30

psmall <- 1.0e-10

pbig <- 1 - 1.0e-7

tmp <- pmin( pbig, pmax( p, psmall ) )

return( log( tmp/(1-tmp) ) )     }
