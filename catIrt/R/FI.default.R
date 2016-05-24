FI.default <-
function( params,                            # parameters over which to calculate
          theta,                             # values/estimates of theta
          type = c("expected", "observed"),  # which information to calculate
          resp = NULL )                      # a response vector/matrix                                
{
	
  return( FI.brm( params, theta, type, resp ) )
    
} # END FI.default FUNCTION