###extract number of locations from Move
if (!isGeneric("n.indiv")) {setGeneric("n.indiv", function(obj) standardGeneric("n.indiv"))}

setMethod("n.indiv", "Move", function(obj){
  return(1)
})

setMethod("n.indiv", ".MoveTrackStack", function(obj){
	  nrow(idData(obj, drop=F))
})

#setMethod("n.locs", ".unUsedRecordsStack", function(obj){
#	  unlist(tapply(obj@trackIdUnUsedRecords, obj@trackIdUnUsedRecords, length))
#})
# 
#setMethod("n.locs", ".unUsedRecords", function(obj){
#	  length( obj@timestampsIdUnUsedRecords)
#})
