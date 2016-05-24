setGeneric("speedSummary", function(x){standardGeneric("speedSummary")})
setMethod("speedSummary", 
	  signature=".MoveTrackSingle",
	  definition=function(x){
		  if(length(seglength(x)>0)){
			  Speed <- speed(x) #meter per sec
			  df  <- data.frame(AverSpeed=mean(Speed, na.rm=TRUE))
			  df$VarSpeed <- var(Speed, na.rm=T)
			  df$MaxSpeed <- max(Speed)
			  return(df)} else {NA}
	  })

setMethod("speedSummary", 
	  signature=".MoveTrackStack", 
	  definition=function(x){
		  lst <- lapply(split(x), speedSummary)
		  df <- do.call("rbind", lst)
		  return(lst)
	  })


setGeneric("speed", function(x){standardGeneric("speed")})
setMethod("speed", 
	  signature=".MoveTrackSingle",
	  definition=function(x){
		  #if(length(seglength(x)>0)){
		  Speed <- (distance(x))/timeLag(x, units="secs") #meter per sec
		  return(Speed)#} else {return(NA)}
	  })

setMethod("speed", 
	  signature=".MoveTrackStack", 
	  definition=function(x){
		  lst <- lapply(split(x), speed)
		  return(lst)
	  })

