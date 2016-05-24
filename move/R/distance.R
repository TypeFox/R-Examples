setGeneric("distanceSummary", function(x){standardGeneric("distanceSummary")})
setMethod("distanceSummary", 
          signature=".MoveTrackSingle",
          definition=function(x){ 
            track <- coordinates(x)
            if(nrow(track)>2){
              Dists <- seglength(x) #vector of all distances in km
              df <- data.frame(TravDist=sum(Dists))    #travel distance in km
              df$MaxDist <- max(Dists) #largest distance
              df$MinDist <- min(Dists) #shortest distance
              if(isLonLat(x)){
                df$FarthDist <- max(spDistsN1(pts=track,pt=track[1,],longlat=T))
              } else {
                df$FarthDist <- max(spDistsN1(pts=track,pt=track[1,],longlat=F))
              } #farthest distance from the start in km
              df$AverDist <- mean(Dists)    #mean distance between relocations in km
              df$SDDist <- sd(Dists)      #standard deviation of distances between relocations
              df$SEDist <- max(as.numeric(spDistsN1(pts=t(as.matrix(track[n.locs(x),])),pt=track[1,],longlat=isLonLat(x)))) #start to end straight distance in km
              return(df)} else {NA}#{warning("Two or less locations.")}
          })

setMethod("distanceSummary", 
          signature=".MoveTrackStack", 
          definition=function(x){
            lst <- lapply(split(x), distanceSummary)
            return(lst)
          })


#setGeneric("distance")#, function(x){standardGeneric("distance")})
setMethod("distance", 
          signature=".MoveTrackSingle",
          definition=function(x){ 
 #           if(n.locs(track)>2){
              Dists <- seglength(x) #vector of all distances in km IF PROJECTION IS LONLAT!!!
#              } else {Dists <- NA}#{warning("Two or less locations.")}
            return(Dists)
          })

setMethod("distance", 
          signature=".MoveTrackStack", 
          definition=function(x){
            lst <- lapply(split(x), distance)
            return(lst)
          })
