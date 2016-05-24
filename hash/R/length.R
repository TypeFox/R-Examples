#  ---------------------------------------------------------------------
# length.R
#   return the number of keys in a hash
#   NB:
#     - This doesn't work: env.profile(x@.xData)$nchains
#  ---------------------------------------------------------------------

setMethod( "length" , "hash" ,
    function(x) 
      length( x@.xData )  
)


