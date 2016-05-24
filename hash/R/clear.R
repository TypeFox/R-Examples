# ---------------------------------------------------------------------
# clear.R
# METHOD: clear
#   clears, rm all key-value pairs from a hash without destroying the 
#   hash
# 
#  TODO:
#   - for large hashes it might be more efficient to re-initialize the 
#     slot than rm the keys on the hash.
#
# ---------------------------------------------------------------------	

setGeneric( "clear", function(x) standardGeneric("clear") )

setMethod( "clear" , "hash" ,
   function(x) rm( list=keys(x), envir=x@.Data )
)


