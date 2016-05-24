setGeneric("trackId", function(x){standardGeneric("trackId")})
setMethod("trackId", 
          signature=".MoveTrackStack",
          definition=function(x){  
		  return(slot(x,'trackId'))
	  })
setMethod("trackId", 
          signature=".unUsedRecordsStack",
          definition=function(x){  
		  return(slot(x,'trackIdUnUsedRecords'))
	  })
