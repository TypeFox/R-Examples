print.SpatialPixelsDataFrame <- function(x, ...)
{
    if (!inherits(x, "SpatialPixelsDataFrame"))
        stop("Non convenient data")
    cat("Object of class \"SpatialPixelsDataFrame\" (package sp):\n\n")
    cat("Grid parameters:\n")
    print(gridparameters(x))
    cat("\nVariables measured:\n")
    print(head(slot(x, "data"), 6))
    if (length(x[[1]])>6)
        cat("...\n")
    cat("\n")
}

print.SpatialPolygonsDataFrame <- function(x, ...)
{
    if (!inherits(x, "SpatialPolygonsDataFrame"))
        stop("Non convenient data")
    cat("Object of class \"SpatialPolygonsDataFrame\" (package sp):\n\n")
    cat("Number of SpatialPolygons: ", length(x[[1]]))
    cat("\n\nVariables measured:\n")
    print(head(slot(x, "data"), 6))
    if (length(x[[1]])>6)
        cat("...\n")
    cat("\n")
}


print.SpatialPixels <- function(x, ...)
{
    if (!inherits(x, "SpatialPixels"))
        stop("Non convenient data")
    cat("Object of class \"SpatialPixels\" (package sp):\n\n")
    cat("Grid parameters:\n")
    print(gridparameters(x))
    cat("\n")
}



print.SpatialGridDataFrame <- function(x, ...)
{
    if (!inherits(x, "SpatialGridDataFrame"))
        stop("Non convenient data")
    cat("Object of class \"SpatialGridDataFrame\" (package sp):\n\n")
    cat("Grid parameters:\n")
    print(gridparameters(x))
    cat("\nVariables measured:\n")
    print(head(slot(x, "data"), 6))
    if (length(x[[1]])>6)
        cat("...\n")
    cat("\n")
}

head.SpatialPoints <- function(x, n=6, ...)
{
    cat("SpatialPoints:\n")
    print(x[1:n,]@coords)
    pst <- paste(strwrap(paste(
                               "Coordinate Reference System (CRS) arguments:",
                               proj4string(x))), collapse="\n")
    cat(pst, "\n")
}



tail.SpatialPoints <- function(x, n=6, ...)
{
    cat("SpatialPoints:\n")
    print(x[(nrow(x)-n):nrow(x),]@coords)
    pst <- paste(strwrap(paste(
                               "Coordinate Reference System (CRS) arguments:",
                               proj4string(x))), collapse="\n")
    cat(pst, "\n")
}



head.SpatialPointsDataFrame <- function(x, n=6, ...)
{
    cc = substring(paste(as.data.frame(t(signif(coordinates(x))))),2,999)
    print(head(data.frame("coordinates" = cc, x@data), n=n, ...))
    pst <- paste(strwrap(paste(
                               "Coordinate Reference System (CRS) arguments:",
                               proj4string(x))), collapse="\n")
    cat(pst, "\n")
}

tail.SpatialPointsDataFrame <- function(x, n=6, ...)
{
    cc = substring(paste(as.data.frame(t(signif(coordinates(x))))),2,999)
    print(tail(data.frame("coordinates" = cc, x@data), n=n, ...))
    pst <- paste(strwrap(paste(
                               "Coordinate Reference System (CRS) arguments:",
                               proj4string(x))), collapse="\n")
    cat(pst, "\n")
}


setMethod("show", "SpatialPixelsDataFrame",
          function(object) print.SpatialPixelsDataFrame(object))

setMethod("show", "SpatialPolygonsDataFrame",
          function(object) print.SpatialPolygonsDataFrame(object))

setMethod("show", "SpatialPixels",
          function(object) print.SpatialPixels(object))

setMethod("show", "SpatialGridDataFrame",
          function(object) print.SpatialGridDataFrame(object))
