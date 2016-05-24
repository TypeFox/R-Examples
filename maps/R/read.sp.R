# transform a SpatialPolygons[DataFrame] into a list of polygons for map()
SpatialPolygons2map <- function(database, namefield=NULL){
  if(!inherits(database,"SpatialPolygons")) stop("database must be a SpatialPolygons[DataFrame] object.")

  region.names <- NULL
  if (inherits(database,"SpatialPolygonsDataFrame") & !is.null(namefield) ) {
    namcol <- lapply(namefield, function(x) which(tolower(names(database)) == tolower(x)))
    if (any(lapply(namcol, length) != 1)) {
      zz <- which(lapply(namcol, length) != 1)
      warning(paste0("database does not (uniquely) contain the field '",namefield[zz],"'."))
    } else {
      zz <- as.data.frame(lapply(database@data[unlist(namcol)], as.character), stringsAsFactors=FALSE)
      region.names <- vapply(1:dim(zz)[1], function(x) paste(zz[x,],collapse=":"),FUN.VALUE="a")
    }
  }
  if (is.null(region.names)) region.names <- unlist(lapply(database@polygons, function(x) x@ID))

  nregions <- length(region.names)

  # count the number of polygons in every "region"
  ngon <- vapply(1:nregions,
                 FUN=function(i) length(database@polygons[[i]]@Polygons),
                 FUN.VALUE=1)
  # if a region contains several polygons, an index is added to the name: "region:n"
  gon.names <- unlist(lapply(1:dim(database)[1], function(i) {
             if (ngon[i]==1) region.names[i]
             else paste(region.names[i],1:ngon[i],sep=":")}))

  # extract all polygon data to a list
  allpoly <- lapply(database@polygons,
                    function(x) lapply(x@Polygons, function(y) y@coords))
## allpoly is a list of lists of Nx2 matrices (not data frames)
## first flatten the list, then add NA to every row, then rbind and remove one NA
#  p1 <- do.call(c, allpoly)
#  p2 <- lapply(p1, function(x) rbind(c(NA,NA),x))
#  p3 <- do.call(rbind,p2)[-1,]
  result <- do.call(rbind, lapply(do.call(c,allpoly),
                                  function(x) rbind(c(NA,NA),x)))[-1,]
  list(x = result[,1], y = result[,2], names = gon.names,
       range = c(range(result[,1], na.rm = TRUE),range(result[,2], na.rm = TRUE)))
}

# transform a SpatialLines[DataFrame] into a list of polylines for map()
SpatialLines2map <- function(database, namefield=NULL){
  if(!inherits(database,"SpatialLines")) stop("database must be a SpatialLines[DataFrame] object.")

  line.names <- NULL
  if (inherits(database,"SpatialLinesDataFrame") & !is.null(namefield) ) {
    namcol <- lapply(namefield, function(x) which(tolower(names(database)) == tolower(x)))
    if (any(lapply(namcol, length) != 1)) {
      zz <- which(lapply(namcol, length) != 1)
      warning(paste0("database does not (uniquely) contain the field '",namefield[zz],"'."))
    } else {
      zz <- as.data.frame(lapply(database@data[unlist(namcol)], as.character), stringsAsFactors=FALSE)
      line.names <- vapply(1:dim(zz)[1], function(x) paste(zz[x,],collapse=":"),FUN.VALUE="a")
    }
  }
  if (is.null(line.names)) line.names <- unlist(lapply(database@lines, function(x) x@ID))

  nlines <- length(line.names)

  # count the number of line segments in every "line"
  nseg <- vapply(1:nlines,
                 FUN=function(i) length(database@lines[[i]]@Lines),
                 FUN.VALUE=1)
  # if a line contains several sub-lines (segments), an index is added to the name: "line:n"
  line.names <- unlist(lapply(1:dim(database)[1], function(i) {
             if (nseg[i]==1) line.names[i]
             else paste(line.names[i],1:nseg[i],sep=":")}))

  # extract all polyline data to a list
  allpoly <- lapply(database@lines,
                    function(x) lapply(x@Lines, function(y) y@coords))
## allpoly is a list of lists of Nx2 matrices (not data frames)
## first flatten the list, then add NA to every row, then rbind and remove one NA
#  p1 <- do.call(c, allpoly)
#  p2 <- lapply(p1, function(x) rbind(c(NA,NA),x))
#  p3 <- do.call(rbind,p2)[-1,]
  result <- do.call(rbind, lapply(do.call(c,allpoly),
                                  function(x) rbind(c(NA,NA),x)))[-1,]
  list(x = result[,1], y = result[,2], names = line.names,
       range = c(range(result[,1], na.rm = TRUE),range(result[,2], na.rm = TRUE)))
}

