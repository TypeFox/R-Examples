KL.default <-
function( params,        # parameters over which to calculate
          theta,         # value of theta
          delta = .1 )   # the indifference region specification
{
  
return( KL.brm( params, theta, delta ) )
 
} # END KL.default FUNCTION