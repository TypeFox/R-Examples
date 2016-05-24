###extract number of locations from Move
if (!isGeneric("n.locs")) {setGeneric("n.locs", function(obj) standardGeneric("n.locs"))}

setMethod("n.locs", "SpatialPointsDataFrame", function(obj){
  return(length(coordinates(obj)[ ,1]))
})

setMethod("n.locs", ".MoveTrackStack", function(obj){
	  unlist(tapply(obj@trackId, obj@trackId, length))
})

#setMethod("n.locs", ".unUsedRecordsStack", function(obj){
#	  unlist(tapply(obj@trackIdUnUsedRecords, obj@trackIdUnUsedRecords, length))
#})
# 
#setMethod("n.locs", ".unUsedRecords", function(obj){
#	  length( obj@timestampsIdUnUsedRecords)
#})
