# Should write a NEW KL function to do integration between theta - delta and theta + delta?
KL <-
function( params,      # parameters over which to calculate
          theta,       # value of theta
          delta = .1 ) # the indifference region specification
{
	  
  UseMethod( "KL" )
  
} # END KL FUNCTION
