setGeneric("summary")
setMethod(f = "summary", 
          signature = c(object = ".UD"), 
          definition = function(object) {
            return(list(Raster_proj=object@crs@projargs, Raster_ext=object@extent, Raster_max_val=maxValue(object), Raster_min_val=minValue(object)))
          })

setMethod("summary", 
          signature=".UDStack", 
          definition=function(object){
            lst <- lapply(split(object), summary)
            return(lst)
          })

setMethod("summary", 
          signature=".MoveTrackSingle", 
          definition=function(object){
            if (!requireNamespace('circular')) 
              stop("You need to install the circular package to proceed") #for angle
            if (!isLonLat(object)) {
              object <- spTransform(object, CRSobj="+proj=longlat")
              warning("The projeciton of the object was changed to longlat inside this function")
            }           
            lst <- list()
            lst <- list(distanceSummary(object), timeSummary(object), speedSummary(object), angleSummary(object))
            names(lst) <- rep(rownames(object@idData), length(lst))
            return(lst)
          })

setMethod("summary", 
          signature=".MoveTrackStack", 
          definition=function(object){
            lapply(split(object), summary)
          })
