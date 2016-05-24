setGeneric("burstId", function(x){standardGeneric("burstId")})
setMethod("burstId", 
          signature=".MoveTrackSingleBurst",
          definition=function(x){  
		  return(slot(x,'burstId'))
	  })
