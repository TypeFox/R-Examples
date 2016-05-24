#  ---------------------------------------------------------------------
#  has-key.R
#  
#  METHOD: has.key( k )  logical
#
#  Returns logical indicating if the hash contains the key  
#  
#  NOTES:
#  - We could implement as a plain old function, but the function then 
#    would have to check based on the signature arguments.  It is just 
#    simpler to use the S4 dispatch mechanism 
# 
#  - See documentation for has-key-methods.Rd and has-key.Rd
#
#  ---------------------------------------------------------------------

setGeneric( 
    "has.key", 
    function( key, hash, ... ) standardGeneric( "has.key" ) 
)


setMethod( 
    "has.key" ,
    signature( "ANY", "hash" ) ,
    function( key, hash, ... ) {
      sapply( key, exists, hash@.Data, inherits=FALSE )
    }
)
