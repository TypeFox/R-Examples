
# this function is internal, not meant for the user

# R code to process user input, then call
# the corresponding function of the C++ object


.sety <- function( z, rWy )  {
  if( NCOL(rWy) > 1 | !is.numeric(rWy) )  
    stop( "'rWy' must be a numeric vector" )
  rWy <- drop(rWy)
  if( !is.null(z$offset) )  rWy <- rWy - z$offset
  rWysorted <- rWy[ order(z$x1) ]
  if( is.vector(z$weights) )  if( any(z$weights==0) )  {
    wsorted <- z$weights[ order(z$x1) ]
    rWysorted <- rWysorted[ wsorted!=0 ]
  }
  (z$CppObj)$sety( rWysorted )
}

