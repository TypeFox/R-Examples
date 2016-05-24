setGeneric("lines") ##need to be different -- 
setMethod("lines", ".MoveTrackSingle", function(x,...){
	  lines(coordinates(x), ...)
})

setMethod("lines", ".MoveTrackStack", function(x,col=NA,...){
	  if(any(is.na(col)))
		  col <- 1:length(unique(x@trackId))
	  if(all(length(col)==length(unique(x@trackId))))
		  col <- col[cumsum(c(1,diff(as.numeric(x@trackId))!=0))] # needs to correspond to points function
	  if(length(col)==1)
		  col<-rep(col, sum(n.locs(x)))
	  if(sum(n.locs(x))!=(length(col)-1) &sum(n.locs(x))!=length(col) )
		  stop('Number of colours does not match either the number of locations/ individuals or is one color')
	  s<-summary(x@trackId)
	  if(any(s==1))
	  {
		  warning('There are/is ',sum(s==1),' individual(s) with only one location for those no line is ploted')
		  s<-s[s!=1]

	  }
	  res <- lapply(names(s), function(Id, x, ...) {
			coords <- coordinates(x)[x@trackId==Id,] 
			graphics::segments(x0=coords[-nrow(coords),1], 
				 y0=coords[-nrow(coords),2], 
				 x1=coords[-1,1], 
				 y1=coords[-1,2], 
				 col=col[x@trackId==Id], ...)
},x=x, ...)
})

setMethod("lines", ".MoveTrackSingleBurst", function(x,col=NA,...){
	  coords <- coordinates(x)


	  if (length(col)==1 && is.na(col)) {
		  col <- x@burstId
	  } else {
		  if (length(col)==length(levels(x@burstId))){
			  col <- col[as.numeric(x@burstId)] # needs to correspond to points function
		  } else {
			  if(length(col)!=1 & length(col)!=n.locs(x))
			  stop("The number of assigned colors is unequal to the number of burst IDs, one or the number of segments")
		  }
	  }
	  if(length(levels(x@burstId))>8) warning("There are more burst IDs than colors (recycling colors).")
	  graphics::segments(x0=coords[-nrow(coords),1], y0=coords[-nrow(coords),2], x1=coords[-1,1], y1=coords[-1,2], col=col, ...)
}) 
