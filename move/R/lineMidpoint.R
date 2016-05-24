setGeneric("lineMidpoint", function(object){standardGeneric("lineMidpoint")})
setMethod(f = "lineMidpoint",
          signature="SpatialPointsDataFrame",
          definition=function(object){
            track <- coordinates(object)
            if(!isLonLat(object)) stop(paste("Coordinates need longlat projection! Projection is:",proj4string(object)))
            if (nrow(track)<2) { ##make the point at the coordinate
              mid <- t(coordinates(track))
            } else {
              if (nrow(track)==2) { ##make the point at the center of the line on a great circle
                mid <- midPoint(coordinates(track)[1,], coordinates(track)[2,])
              }
              if (nrow(track)>2){
                dists <- seglength(object)
                cumsum <- cumsum(dists)
                min <- which(cumsum(dists)>0.5*sum(dists))[1]
                if (sum(dists)==0) {min <- 0} else {min <- which(cumsum>0.5*sum(dists))[1]}
                if (min==0) {mid <- coordinates(track)[1,]} else {
                 if (min==1){
                    prop <- ((0.5*sum(dists)))/(cumsum(dists)[min])
                    #####ADDING CALCULATION ALONG GREAT CIRCLE USING DESTPOINT FROM GEOSPHERE
                    mid <- destPoint(p=coordinates(track)[min,], 
                                     b=bearing(coordinates(track)[min,], coordinates(track)[min+1,]), 
                                     d=prop*dists[min])
                  }else{        
                    prop <- (sum(dists)/2-(cumsum(dists)[min-1]))/(cumsum(dists)[min]-(cumsum(dists)[min-1]))
                    #####ADDING CALCULATION ALONG GREAT CIRCLE USING DESTPOINT FROM GEOSPHERE
                    mid <- destPoint(p=coordinates(track)[min,], 
                                     b=bearing(coordinates(track)[min,], coordinates(track)[min+1,]), 
                                     d=prop*dists[min]) #seglength returns km,for lonlat=T,destPoint uses m
                  }
                }
              }
            }
            midSP <- SpatialPoints(coords=as.data.frame(x=mid, row.names=1:nrow(mid)), proj4string=CRS(proj4string(object)))
          })
