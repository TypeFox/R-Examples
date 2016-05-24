# ------------------------------------------------------------------------------
# Methods for class 'SegSpatial'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Coercion methods
# ------------------------------------------------------------------------------
setAs("SegSpatial", "SpatialPoints", 
      function(from) {
        validObject(from)
        SpatialPoints(coords = from@coords, proj4string = from@proj4string)
      })

setAs("SegSpatial", "SpatialPointsDataFrame", 
      function(from) {
        validObject(from)
        SpatialPointsDataFrame(coords = from@coords, 
                               data = data.frame(from@data),
                               proj4string = from@proj4string)
      })

setAs("SegSpatial", "SpatialPixelsDataFrame", 
      function(from) {
        validObject(from)
        SpatialPixelsDataFrame(points = from@coords, 
                               data = data.frame(from@data),
                               proj4string = from@proj4string)
      })

#setAs("SegSpatial", "SegLocal", 
#      function(from) {
#        validObject(from)
#        SegLocal(coords = from@coords, data = from@data, env = from@env, 
#                 proj4string = from@proj4string)
#      })

setAs("SegSpatial", "list", 
      function(from) {
        validObject(from)
        list(d = from@d, r = from@r, h = from@h, p = from@p)
      })

as.list.SegSpatial <- function(x, ...) {
  validObject(x)
  list(d = x@d, r = x@r, h = x@h, p = x@p)
}

# ------------------------------------------------------------------------------
# Printing
# ------------------------------------------------------------------------------
setMethod("show", signature(object = "SegSpatial"), function(object) {
  validObject(object)
  cat("\n\tReardon and O'Sullivan's spatial segregation measures\n\n")
  
  cat("Dissimilarity (D)     : ")
  if (length(object@d) > 0)
    cat(round(object@d, 4), "\n")
  else
    cat("-\n")
  
  cat("Relative diversity (R): ")
  if (length(object@r) > 0)
    cat(round(object@r, 4), "\n")
  else
    cat("-\n")
  
  cat("Information theory (H): ")            
  if (length(object@h) > 0)
    cat(round(object@h, 4), "\n")
  else
    cat("-\n")
  
  cat("Exposure/Isolation (P): ")
  if (length(object@p) > 0) {
    cat("\n")
    if (is.null(colnames(object@p)))
      colnames(object@p) <- paste("Group", 1:ncol(object@p))
    if (is.null(rownames(object@p)))
      rownames(object@p) <- paste("Group", 1:nrow(object@p))
    print(object@p)
    cat("--\n")
    cat("The exposure/isolation matrix should be read horizontally.\n")
    cat("Read 'help(spseg)' for more details.\n")
  } else {
    cat("-\n")
  }
})

print.SegSpatial <- function(x, digits = getOption("digits"), ...) {
  validObject(x)
  cat("\n\tReardon and O'Sullivan's spatial segregation measures\n\n")
  
  cat("Dissimilarity (D)     : ")
  if (length(x@d) > 0)
    cat(round(x@d, digits), "\n")
  else
    cat("-\n")
  
  cat("Relative diversity (R): ")
  if (length(x@r) > 0)
    cat(round(x@r, digits), "\n")
  else
    cat("-\n")
  
  cat("Information theory (H): ")            
  if (length(x@h) > 0)
    cat(round(x@h, digits), "\n")
  else
    cat("-\n")
  
  cat("Exposure/Isolation (P): ")
  if (length(x@p) > 0) {
    cat("\n")
    if (is.null(colnames(x@p)))
      colnames(x@p) <- paste("Group", 1:ncol(x@p))
    if (is.null(rownames(x@p)))
      rownames(x@p) <- paste("Group", 1:nrow(x@p))
    print(x@p, digits, ...)
    cat("--\n")
    cat("The exposure/isolation matrix should be read horizontally.\n")
    cat("Read 'help(spseg)' for more details.\n")
  } else {
    cat("-\n")
  }
}

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
setMethod("spplot", signature(obj = "SegSpatial"), function(obj, ...) {
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

# ------------------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------------------
