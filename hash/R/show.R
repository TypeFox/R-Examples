# -----------------------------------------------------------------------------
# show.r
# 
# METHODS: show 
#  The default method on the class.  Perhaps this should return the 
#  length.
#
#  See Also: print
# -----------------------------------------------------------------------------
setMethod( "show" , "hash" ,
	function(object) cat(format(object))
)
    
    
    
