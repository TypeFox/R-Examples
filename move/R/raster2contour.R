setGeneric("raster2contour", function(x, ...) {standardGeneric("raster2contour")})  
setMethod(f = "raster2contour", 
          signature = c(x = ".UD"),
          definition = function(x, ...) {
            vol<-getVolumeUD(x)
            rasterToContour(vol, ...)
          })

setMethod(f = "raster2contour", 
          signature = c(x = ".UDStack"),
          definition = function(x, ...) {
            xx<-split(x)
            tmp<-lapply(xx, raster2contour, ...=...)
      	    x<-mapply(spChFIDs, tmp, mapply(paste,lapply(tmp,row.names), names(tmp), SIMPLIFY=F))
      	    tmp<-do.call('rbind',mapply('[[<-',MoreArgs=list(i='individual.local.identifier'),x=x, value=names(tmp)))
      	    return(tmp)
          })
