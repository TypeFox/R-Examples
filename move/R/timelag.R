###extract time.lag from Move
#if (!isGeneric("time.lag")) {
#	setGeneric("time.lag", function(x, ...) standardGeneric("time.lag"))
#}
#setMethod("time.lag", 
#	  ".MoveTrackSingle", 
#	  function(x, ...){
#	          .Deprecated('timeLag', package='move', msg='This function will be removed in upcoming versions of the move package, please use timeLag')
#		  return(as.numeric(diff(timestamps(x)),...)) 
#	  }
#	  )
#
#setMethod("time.lag", ".MoveTrackStack", function(x, units,...){
#		  .Deprecated('timeLag', package='move', msg='This function will be removed in upcoming versions of the move package, please use timeLag')
#		  if(missing(units)){
#			  warning('Units not specified this could lead to different units for the time differences between individuals')
#			  return(lapply(lapply(split(timestamps(x), x@trackId), diff),as.numeric,...))
#		  }else{
#			  return(lapply(lapply(split(timestamps(x), x@trackId), diff),as.numeric,units=units,...))
#		  }
#	  })

###extract time.lag from Move
if (!isGeneric("timeLag")) {
	setGeneric("timeLag", function(x, ...) standardGeneric("timeLag"))
}
setMethod("timeLag", 
	  ".MoveTrackSingle", 
	  function(x, ...){
		  return(as.numeric(diff(timestamps(x)),...)) 
	  }
	  )

setMethod("timeLag", ".MoveTrackStack", function(x, units,...){
		  if(missing(units)){
			  warning('Units not specified this could lead to different units for the time differences between individuals')
			  return(lapply(lapply(split(timestamps(x), x@trackId), diff),as.numeric,...))
		  }else{
			  return(lapply(lapply(split(timestamps(x), x@trackId), diff),as.numeric,units=units,...))
		  }
	  })
