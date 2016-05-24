if (!isGeneric("spCbind"))
	setGeneric("spCbind", function(obj, x)
		standardGeneric("spCbind"))

cbindSpatialPointsDataFrame <- function(obj, x) {
    x0 <- slot(obj, "data")
    if (nrow(x0) != nrow(x)) stop("different numbers of rows")
    cx <- data.frame(x0, x)
    SpatialPointsDataFrame(as(obj, "SpatialPoints"), data=cx)
}

cbindSpatialPointsDataFramev <- function(obj, x) {
    x0 <- slot(obj, "data")
    if (nrow(x0) != length(x)) stop("different numbers of rows")
    nx <- deparse(substitute(x))
    x <- as.data.frame(x)
    names(x) <- nx
    cx <- data.frame(x0, x)
    SpatialPointsDataFrame(as(obj, "SpatialPoints"), data=cx)
}

setMethod("spCbind", signature(obj="SpatialPointsDataFrame", x="data.frame"), 
    cbindSpatialPointsDataFrame)

setMethod("spCbind", signature(obj="SpatialPointsDataFrame", x="vector"), 
    cbindSpatialPointsDataFramev)

cbindSpatialLinesDataFrame <- function(obj, x) {
    x0 <- slot(obj, "data")
    if (nrow(x0) != nrow(x)) stop("different numbers of rows")
    if (!isTRUE(all.equal(row.names(x0), row.names(x))))
        stop("row names not identical")
    cx <- data.frame(x0, x)
    SpatialLinesDataFrame(as(obj, "SpatialLines"), data=cx)
}

cbindSpatialLinesDataFramev <- function(obj, x) {
    x0 <- slot(obj, "data")
    if (nrow(x0) != length(x)) stop("different numbers of rows")
    nx <- deparse(substitute(x))
    x <- as.data.frame(x)
    names(x) <- nx
    cx <- data.frame(x0, x)
    SpatialLinesDataFrame(as(obj, "SpatialLines"), data=cx)
}

setMethod("spCbind", signature(obj="SpatialLinesDataFrame", x="data.frame"), 
    cbindSpatialLinesDataFrame)

setMethod("spCbind", signature(obj="SpatialLinesDataFrame", x="vector"), 
    cbindSpatialLinesDataFramev)

cbindSpatialPolygonsDataFrame <- function(obj, x) {
    x0 <- slot(obj, "data")
    if (nrow(x0) != nrow(x)) stop("different numbers of rows")
    if (!isTRUE(all.equal(row.names(x0), row.names(x))))
        stop("row names not identical")
    cx <- data.frame(x0, x)
    SpatialPolygonsDataFrame(as(obj, "SpatialPolygons"), data=cx)
}

cbindSpatialPolygonsDataFramev <- function(obj, x) {
    x0 <- slot(obj, "data")
    if (nrow(x0) != length(x)) stop("different numbers of rows")
    nx <- deparse(substitute(x))
    x <- as.data.frame(x)
    names(x) <- nx
    cx <- data.frame(x0, x)
    SpatialPolygonsDataFrame(as(obj, "SpatialPolygons"), data=cx)
}

setMethod("spCbind", signature(obj="SpatialPolygonsDataFrame", x="data.frame"), 
    cbindSpatialPolygonsDataFrame)

setMethod("spCbind", signature(obj="SpatialPolygonsDataFrame", x="vector"), 
    cbindSpatialPolygonsDataFramev)


if (!isGeneric("spRbind"))
	setGeneric("spRbind", function(obj, x)
		standardGeneric("spRbind"))

rbindSpatialPoints <- function(obj, x) {
    if (!isTRUE(all.equal(proj4string(obj), proj4string(x))))
        stop("coordinate reference systems differ")
    crds <- rbind(coordinates(obj), coordinates(x))
    SpatialPoints(crds, proj4string=CRS(proj4string(obj)))
}

setMethod("spRbind", signature(obj="SpatialPoints", x="SpatialPoints"),
    rbindSpatialPoints)

rbindSpatialPointsDataFrame <- function(obj, x) {
    SP <- spRbind(as(obj, "SpatialPoints"), as(x, "SpatialPoints"))
#    df <- rbind(as(obj, "data.frame"), as(x, "data.frame"))
#    stopped adding coordinates as variables; Steve Eick 100117
    df <- rbind(slot(obj, "data"), slot(x, "data"))
    SpatialPointsDataFrame(SP, data=df)
}

setMethod("spRbind", signature(obj="SpatialPointsDataFrame", 
    x="SpatialPointsDataFrame"), rbindSpatialPointsDataFrame)

rbindSpatialLines <- function(obj, x) {
    if (!isTRUE(all.equal(proj4string(obj), proj4string(x))))
        stop("coordinate reference systems differ")
    ll1 <- slot(obj, "lines")
    ll2 <- slot(x, "lines")
    ID1 <- sapply(ll1, function(x) slot(x, "ID"))
    ID2 <- sapply(ll2, function(x) slot(x, "ID"))
    if (length(c(ID1, ID2)) > length(unique(c(ID1, ID2))))
        stop("non-unique line IDs")
    LL <- c(ll1, ll2)
    SpatialLines(LL, proj4string=CRS(proj4string(obj)))
}

setMethod("spRbind", signature(obj="SpatialLines", x="SpatialLines"),
    rbindSpatialLines)

rbindSpatialLinesDataFrame <- function(obj, x) {
    SL <- spRbind(as(obj, "SpatialLines"), as(x, "SpatialLines"))
    df <- rbind(as(obj, "data.frame"), as(x, "data.frame"))
    SpatialLinesDataFrame(SL, data=df)
}

setMethod("spRbind", signature(obj="SpatialLinesDataFrame", 
    x="SpatialLinesDataFrame"), rbindSpatialLinesDataFrame)


rbindSpatialPolygons <- function(obj, x) {
    if (!isTRUE(all.equal(proj4string(obj), proj4string(x))))
        stop("coordinate reference systems differ")
    pl1 <- slot(obj, "polygons")
    pl2 <- slot(x, "polygons")
    ID1 <- sapply(pl1, function(x) slot(x, "ID"))
    ID2 <- sapply(pl2, function(x) slot(x, "ID"))
    if (length(c(ID1, ID2)) > length(unique(c(ID1, ID2))))
        stop("non-unique polygon IDs")
    PL <- c(pl1, pl2)
    SpatialPolygons(PL, proj4string=CRS(proj4string(obj)))
}

setMethod("spRbind", signature(obj="SpatialPolygons", x="SpatialPolygons"),
    rbindSpatialPolygons)

rbindSpatialPolygonsDataFrame <- function(obj, x) {
    SP <- spRbind(as(obj, "SpatialPolygons"), as(x, "SpatialPolygons"))
    df <- rbind(as(obj, "data.frame"), as(x, "data.frame"))
    SpatialPolygonsDataFrame(SP, data=df)
}

setMethod("spRbind", signature(obj="SpatialPolygonsDataFrame", 
    x="SpatialPolygonsDataFrame"), rbindSpatialPolygonsDataFrame)





