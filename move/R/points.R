setGeneric("points")
setMethod("points", ".MoveTrackStack", function(x,col=NA,...){
  if(any(is.na(col)))
    col <- 1:length(unique(x@trackId))
	  if(all(length(col)==length(unique(x@trackId))))
		  col <- col[cumsum(c(1,diff(as.numeric(x@trackId))!=0))]# needs to correspond to lines function
  if(sum(n.locs(x))!=length(col) & length(col)!=1)
	  stop('Number of colours does not match either the number of locations/ individuals or is one color')
  points(coordinates(x), col=col,...)
})

setMethod("points", ".MoveTrackSingle", function(x,...){
  points(coordinates(x), ...)
})
