# -----------------------------------------------------------------------------
# CONSTRUCTOR: hash
#   Takes an optional 1 or two parameter
#   DEPENDS on method set
# -----------------------------------------------------------------------------
hash <- function( ... ) {

  li <- list(...)  

  # INITIALIZE A NEW HASH   
  h <- new( 
    "hash" , 
     new.env( hash = TRUE , parent=emptyenv() )  
  )

  if ( length(li) >  0  ) { 
    if( length(li) > 0 ) .set( h, ... )
  }

  return(h)

}
	


