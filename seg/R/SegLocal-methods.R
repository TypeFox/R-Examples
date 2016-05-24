# ------------------------------------------------------------------------------
# Methods for class 'SegLocal'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Coercion methods
# ------------------------------------------------------------------------------
setAs("SegLocal", "SpatialPoints", 
      function(from) {
        validObject(from)
        SpatialPoints(coords = from@coords, proj4string = from@proj4string)
      })

setAs("SegLocal", "SpatialPointsDataFrame", 
      function(from) {
        validObject(from)
        SpatialPointsDataFrame(coords = from@coords, 
                               data = data.frame(from@env),
                               proj4string = from@proj4string)
      })

setAs("SegLocal", "SpatialPixelsDataFrame", 
      function(from) {
        validObject(from)
        SpatialPixelsDataFrame(points = from@coords, 
                               data = data.frame(from@env),
                               proj4string = from@proj4string)
      })

setAs("SpatialPointsDataFrame", "SegLocal", 
      function(from) {
        SegLocal(coords = from@coords, data = from@data, env = from@data, 
                 proj4string = from@proj4string)
      })

setAs("SpatialPolygonsDataFrame", "SegLocal", 
      function(from) {
        SegLocal(coords = coordinates(from), data = from@data, 
                 env = from@data, proj4string = from@proj4string)
      })

# ------------------------------------------------------------------------------
# Printing
# ------------------------------------------------------------------------------
setMethod("show", signature(object = "SegLocal"), function(object) {
  validObject(object)
  cat("Class                 :", class(object), "\n")
  cat("Number of data points :", nrow(object@coords), "\n")
  cat("Number of data columns:", ncol(object@data), "\n")
  cat("Projection            :", object@proj4string@projargs, "\n")
  cat("Slot names            :", slotNames(object), "\n")
})

print.SegLocal <- function(x, ...) {
  validObject(x)
  cat("Class                 :", class(x), "\n")
  cat("Number of data points :", nrow(x@coords), "\n")
  cat("Number of data columns:", ncol(x@data), "\n")
  cat("Projection            :", x@proj4string@projargs, "\n")
  cat("Slot names            :", slotNames(x), "\n")
}

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
plot.SegLocal <- function(x, which.col, main, ...) {
  validObject(x)
  xx <- x@coords[,1]
  yy <- x@coords[,2]
  env <- x@env
  
  if (missing(which.col))
    which.col <- 1:ncol(env)
  numPlot <- length(which.col) 
  if (missing(main))
    main <- paste("Data", 1:numPlot)
  else if (length(main) < numPlot)
    main <- rep(main, ceiling(numPlot/length(main)))
  
  for (i in 1:numPlot) {
    qq <- quantile(env[,i], probs = c(0.2, 0.4, 0.6, 0.8))
    brks <- length(qq) + 1
    size <- rep(brks, nrow(env))
    for (j in 1:brks)
      size[which(env[,i] <= qq[brks - j])] <- brks - j
    plot(x = xx, y = yy, cex = size, main = main[i], ...)
  }
}

points.SegLocal <- function(x, which.col = 1, ...) {
  validObject(x)
  xx <- x@coords[,1]
  yy <- x@coords[,2]
  env <- x@env
  
  if (length(which.col) > 1)
    warning("'which.col' has a length of > 1", call. = FALSE)
  i <- which.col[1] 
  
  qq <- quantile(env[,i], probs = c(0.2, 0.4, 0.6, 0.8))
  brks <- length(qq) + 1
  size <- rep(brks, nrow(env))
  for (j in 1:brks)
    size[which(env[,i] <= qq[brks - j])] <- brks - j
  points(x = xx, y = yy, cex = size, ...)
}

setMethod("spplot", signature(obj = "SegLocal"), function(obj, ...) {
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
summary.SegLocal <- function(object, ...) {
  validObject(object)
  msg <- paste("An object of class \"", class(object), "\"\n", sep = "")
  cat(msg)
  cat("Coordinates:\n")
  tmp <- t(apply(object@coords, 2, range))
  rownames(tmp) <- c("x", "y")
  colnames(tmp) <- c("min", "max")
  print(tmp)
  if (is.na(object@proj4string@projargs)) {
    cat("Is projected: FALSE\n")
  } else {
    cat("Is projected: TRUE\n")
    cat("proj4string : [", object@proj4string@projargs, "]\n", sep = "")
  }
    
  cat("\nData values (%):\n")
  tmp <- t(apply(object@data, 1, function(z) z/sum(z))) * 100
  tmp <- apply(tmp, 2, function(z) summary(z, ...))
  colnames(tmp) <- colnames(object@data)
  print(tmp)

  cat("\nLocal environment composition (%):\n")
  tmp <- t(apply(object@env, 1, function(z) z/sum(z))) * 100
  tmp <- apply(tmp, 2, function(z) summary(z, ...))
  colnames(tmp) <- colnames(object@env)
  print(tmp)
}

update.SegLocal <- function(object, coords, data, env, proj4string, ...) {
  validObject(object)
  if (missing(coords))
    coords <- object@coords  
  if (missing(data))
    data <- object@data
  if (missing(env))
    env <- object@env
  if (missing(proj4string))
    proj4string <- object@proj4string
              
  SegLocal(coords, data, env, proj4string)
}

# Methods that are not so useful (removed on 23 December 2013) ...
#
# as.list.SegLocal <- function(x, ...) {
#   validObject(x)
#   list(coords = x@coords, data = x@data, env = x@env, 
#        proj4string = x@proj4string)
# }
# 
# setAs("list", "SegLocal",
#       function(from) {
#         if (is.null(from$proj4string))
#           SegLocal(coords = from$coords, data = from$data, env = from$env)
#         else
#           SegLocal(coords = from$coords, data = from$data, env = from$env, 
#                    proj4string = from$proj4string)
#       })
# setAs("SegLocal", "list", 
#       function(from) {
#         validObject(from)
#         list(coords = from@coords, data = from@data, env = from@env, 
#              proj4string = from@proj4string)
#       })
#
# "[[.SegLocal" <- function(i, ...) {
#   validObject(x)
#   slotnames <- slotNames(x)
# 
#   if (is.numeric(i)) {
#     i <- as.integer(i)
#     if (i > length(slotnames))
#       chosen <- NULL
#     else {
#       chosen <- slotnames[i]
#       chosen <- paste("x@", chosen, sep = "")
#       chosen <- eval(parse(text = chosen))
#     }
#   }
#   
#   else if (is.character(i)) {
#     chosen <- paste("x@", i, sep = "")
#     chosen <- eval(parse(text = chosen))
#   }
#   
#   else {
#     chosen <- NULL
#   }
# }
