
setGeneric("turnAngleGc", function(x){standardGeneric("turnAngleGc")})
setMethod("turnAngleGc", 
          signature=".MoveTrackSingle",
          definition=function(x){
		  if(!isLonLat(x))
			  warning('turnAngleGc is probably not a valid calculation on this projection')
		 fb<- finalBearing(x[-n.locs(x),], x[-1,])[-(n.locs(x)-1)]
		 b<- bearing(x[-n.locs(x),], x[-1,])[-1]
		 t<-(b-((fb)))
		 return(((t+180)%%360)-180)
          })
setMethod("turnAngleGc", 
          signature=".MoveTrackStack",
          definition=function(x){
            if(!isLonLat(x))
              warning('turnAngleGc is probably not a valid calculation on this projection')
            return(lapply(split(x), turnAngleGc))
          })