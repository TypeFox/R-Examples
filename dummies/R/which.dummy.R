# -----------------------------------------------------------------------------
# FUNCTIONS: which.dummy which dummies
#   Which variables are dummy variables
#  
#  TODO: 
#   - allow for multiple names.
# -----------------------------------------------------------------------------

which.dummy <- function(data, name=NULL) {

  indexes <- integer()  

  if( ! is.null(name) ) {
     indexes <- attr( data, 'dummies' )[[name]] 
  } else {  

    if( is.null( attr( data, 'dummies' ) ) )
      stop( "Data does not appear to have dummy variables." )

    for( name in names( attr( data, 'dummies' ) ) )
      indexes <- append( indexes, attr( data, 'dummies')[[name]] )

      # indexes <- sapply( attr( data, 'dummies' ), I, USE.NAMES=F ) 

  }

  return( sort( as.integer( indexes ) ) )
}

