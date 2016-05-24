# ---------------------------------------------------------------------
# CLASS: hash
# 
# n.b.: The name of this class from the Open Data R Style Guide which 
# would have the CLASS named 'Hash'.  We use hash since the goal is to 
# emulate a native class that is missing from the R Specification.
# 
# TODO:
#
# CONTAINS: hash class and hash accessors ($ [ [[)
#
# --------------------------------------------------------------------- 
setClass( 
  'hash', 
  contains = 'environment' ,
)


# ----------------------- ACCESSOR METHODS ------------------------------

  # --------------------------------------------------------------------- 
  # METHOD: [ (hash slice)
  #   The [ method provides for the subseting of the object and 
  #   extracting a copy of the slice.
  #
  #   Notes:
  #    - Uses 'mget' internally for speed.  Provides access to the hash.  
  #    - We do not use the .set method for speed.   
  # --------------------------------------------------------------------- 

  setMethod( 
       '[' , 
        signature( x="hash", i="ANY", j="missing", drop = "missing") ,  
        function( 
          x,i,j, ... , 
          # na.action = 
          #  if( is.null( getOption('hash.na.action') ) ) NULL else 
          #  getOption('hash.na.action') , 
          drop 
        ) {
  
          .h <- hash() # Will be the new hash
          for( k in i ) assign( k, get(k,x), .h@.Data )
            
          return(.h)

        }
  )

#  system.time( for( i in 1:10 ) for( ke in kes ) ha[ ke ]  )

  # NB. A slice without any arguments, by definition returns the hash itself
  setMethod( '[', signature( 'hash', 'missing', 'missing', 'missing' ),
    function(x,i,j, ..., drop ) {
      return( x )                  
    }
  )



# --------------------------------------------------------------------- 
# METHOD: [<-, Hash Slice Replacement Method
# WHAT DO WE DO IF WE HAVE A DIFFERENT NUMBER OF KEYS AND VALUES?
#   This should implement a hash slice.
#   NB.  Although we would like to use assign directly, we use set 
#        because it deals with the ambiguity of the lengths of the 
#        key and value vectors.
# --------------------------------------------------------------------- 

setReplaceMethod( '[', c(x ="hash", i="ANY" ,j="missing", value="ANY") ,
	function( x, i, ...,  value ) {
	  .set( x, i, value, ...  )  
	  return( x )
    }
)




# hash[ key ] <- NULL : Removes key-value pair from hash
setReplaceMethod( '[', c(x="hash", i="ANY", j="missing", value="NULL") ,
    function( x, i, ...,  value ) {
      del( i, x )
      return( x )
    }
)
  


# TEST:
# h[ "foo" ] <- letters # Assigns letters, a vector to "foo"
# h[ letters ] <- 1:26
# h[ keys ] <- value
# h[ 'a' ] <- NULL 


# ---------------------------------------------------------------------
# $ -- DEPRECATED
#   This is deprecated since '$' is defined on environments and 
#   environments can be inherited in objects
#
# ---------------------------------------------------------------------

# SPECIAL CASE: NULL value
#   When assign a null values to a hash the key is deleted. It is 
#   idiomatic when setting a value to NULL in R that that value is
#   removed from a list or environment. 
#   
#   If R's behavior changes this will go away.
#   It is interesting to note that [[ behaves this way
#
setReplaceMethod( '$', c( x="hash", value="NULL"),
  function(x, name, value) {
    remove( list=name, envir=x@.xData )
    x
  }
)


# ---------------------------------------------------------------------
# [[ -- DEPRECATED:
#   This is deprecated since this is handled by R natively.
#   Return single value, key,i, is a name/gets interpretted.
# 
#   NB: We no longer use .get.
# ---------------------------------------------------------------------

setReplaceMethod( '[[', c(x="hash", i="ANY", j="missing", value="ANY") ,
  function(x,i,value) {
    assign( i, value, x@.xData )
    return( x )
  }
)


# CASE: hash$value <- NULL
#   Deletes the value  
setReplaceMethod( '[[', c(x="hash", i="ANY", j="missing", value="NULL") ,
  function(x,i,value) {
    rm( list=i, envir=x@.xData )
    return( x )
  }
)


# ---------------------------------------------------------------------
# MISC. FUNCTIONS
# ---------------------------------------------------------------------

is.hash <- function(x) is( x, "hash" )

as.list.hash <- function(x, all.names=FALSE, ...) 
  as.list( x@.Data, all.names, ... )

is.empty <- function(x) { 
    if( class(x) != 'hash' ) stop( "is.empty only works on hash objects" )
    if( length(x) == 0 ) TRUE else FALSE  
}

