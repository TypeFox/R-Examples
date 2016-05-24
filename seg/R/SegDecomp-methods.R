# ------------------------------------------------------------------------------
# Methods for class 'SegDecomp'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Coercion methods
# ------------------------------------------------------------------------------
setAs("SegDecomp", "SpatialPoints", 
      function(from) {
        validObject(from)
        SpatialPoints(coords = from@coords, proj4string = from@proj4string)       
      })

setAs("SegDecomp", "SpatialPointsDataFrame", 
      function(from) {
        validObject(from)
        SpatialPointsDataFrame(coords = from@coords, 
                               data = data.frame(from@data),
                               proj4string = from@proj4string)      
      })

setAs("SegDecomp", "SpatialPixelsDataFrame", 
      function(from) {
        validObject(from)
        SpatialPixelsDataFrame(points = from@coords, 
                               data = data.frame(from@data),
                               proj4string = from@proj4string)     
      })

setAs("SegDecomp", "vector", 
      function(from) {
        validObject(from)
        as.vector(from@d)
      })

as.vector.SegDecomp <- function(x, ...) {
  validObject(x)
  as.vector(x@d)
}

# ------------------------------------------------------------------------------
# Printing
# ------------------------------------------------------------------------------
setMethod("show", signature(object = "SegDecomp"), function(object) {
  validObject(object)
  print(sum(object@d))
  cat("\nSubcomponents:\n")
  tmp <- t(data.frame(object@d))
  rownames(tmp) <- ""
  colnames(tmp) <- c("Location", "Composition", "Qualitative")
  print(tmp)
})

print.SegDecomp <- function(x, digits = getOption("digits"), ...) {
  validObject(x)
  print(sum(x@d), digits = digits)
  cat("\nSubcomponents:\n")
  tmp <- t(data.frame(x@d))
  rownames(tmp) <- ""
  colnames(tmp) <- c("Location", "Composition", "Qualitative")
  print(tmp, digits = digits)
}

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
setMethod("spplot", signature(obj = "SegDecomp"), function(obj, ...) {
  validObject(obj)
  spO <- try(as(obj, "SpatialPixelsDataFrame"), silent = TRUE)
  if (class(spO) == "try-error") {
    warning("failed to convert 'obj' to SpatialPixelsDataFrame", call. = FALSE)
    
    spO <- try(as(obj, "SpatialPointsDataFrame"), silent = TRUE)
    if (class(spO) == "try-error")
      stop("failed to convert 'obj' to SpatialPointsDataFrame", call. = FALSE)
  }
  
  spplot(spO, ...)
})
