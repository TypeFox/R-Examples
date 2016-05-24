
coords2json <- function(x, type, holes=NULL, order=NULL){   #x=coords; holes=booleans vector; type=class
  if(type=="splpoints"||type=="splicons"){
    res <- paste("[", x[,1], ", ", x[,2], "]", sep="")
  }
  if(type=="spllines"){
    res <- lapply(x, function(x) lapply(x, proclines))
    res <- lapply(res, function(x) do.call("paste", c(x, sep=", ")))
    res <- lapply(res, function(x) paste("[", x, "]", sep=""))
  }
  if(type=="splpolygons"){
    if(is.null(holes)){
      res <- lapply(x, function(x) lapply(x, procpolys))
      res <- lapply(res, function(x) do.call("paste", c(x, sep=", ")))
      res <- lapply(res, function(x) paste("[", x, "]", sep=""))
    } else {
      res <- lapply(x, function(x) lapply(x, procpolys))
      res <- mapply(function(a, b) a[b], a=res, b=order)
      holes <- mapply(function(a, b) a[b], a=holes, b=order)
      res <- mapply(holify, a=holes, b=res)
      res <- lapply(res, function(x) do.call("paste", c(x, sep=", ")))
      res <- lapply(res, function(x) paste("[", x, "]", sep=""))
    }
  }
  return(res)
}

holify <- function(a, b){
  for(i in 1:length(a)){
    if (a[i]==T && i>1){
      b[i-1] <- substr(b[i-1], 1, nchar(b[i-1])-1)
      b[i] <- substr(b[i], 2, nchar(b[i]))
    }
  }
  return(b)
}

polycoords <- function(x){  #Extrait les coordonnees d'un objet SpatialPolygons
  res <- x@polygons
  res <- lapply(res, function(x) x@Polygons)
  res <- lapply(res, function(x) lapply(x, coordinates))
  return(res)
}

polyholes <- function(x){  #Extrait les slots 'hole' d'un objet SpatialPolygons
  res <- x@polygons
  res <- lapply(res, function(x) x@Polygons)
  res <- lapply(res, function(x) lapply(x, function(x) x@hole))
  return(res)
}

polyorder <- function(x){  #Extrait les slots 'plotOrder' d'un objet SpatialPolygons
  res <- x@polygons
  res <- lapply(res, function(x) x@plotOrder)
  return(res)
}

proclines <- function(x) {
  res <- paste("[", x[,1], ", ",x[,2], "]", sep="", collapse=", ")
  res <- paste("[", res, "]", sep="")
}

procpolys <- function(x) {
  res <- paste("[", x[,1], ", ",x[,2], "]", sep="", collapse=", ")
  res <- paste("[[", res, "]]", sep="")
}

pngasp <- function(x){ #x is an sp bounding box.
  xlim <- x[1,]
  ylim <- x[2,]
  asp <- (diff(ylim)/diff(xlim))/cos((mean(ylim) * pi)/180)
  names(asp) <- NULL
  return(asp)
}


class2var <- function(x){
  switch(x,
         "basemap"="BaseMap",
         "splpoints"="Points",
         "splicons"="Icons",
         "spllines"="Lines",
         "splpolygons"="Polygons",
         "splgrid"="Raster")
}

xvarnames <- function(x){ #x est une liste d'objets spl
  xclass <- sapply(x, function(x) return(class(x)))
  xname <- sapply(x, function(x) return(x$name))
  xvarclass <- sapply(xclass, class2var)
  xvarname <- paste(safeVar(xname), xvarclass, sep="")
  res <- data.frame(xclass, xname, xvarname)
  return (res)
}

safeVar <- function(x){
  x <- gsub("[^A-Za-z0-9]", "", x)
  test1 <- grepl("[0-9]", substr(x, 1, 1))
  substr(x[test1], 1, 1) <- letters[as.numeric(substr(x[test1], 1, 1))+1]

  if(any(x == "")){
    stop("Incorrect value for name")
  }
  return(x)
}

cleanDepsub <- function(x){
  x <- paste(x, collapse="")
  x <- substr(x, 5, nchar(x))
  x <- gsub("\\(", "", x)
  x <- gsub("\\)", "", x)
  x <- gsub(" ", "", x)
  x <- strsplit(x, ",")
  return(x)
}

#' Convert colors to hexadecimal format
#'
#' This function converts any color of the \R system to hexadecimal format.
#' @param x a vector of any of the three kinds of \R color specifications.
#' @param alpha.channel logical. Sould an alpha channel included in the output?  Default is \code{FALSE}.
#' @param alpha a vector of numeric values in \eqn{[0, 1]}. Recycled if necessary.
#' @param charstring logical. Sould resulting elements be surrounded by quotation marks? Default is \code{TRUE}.
#' 
#' @return A character vector of hexadecimal values.
#' 
#' @seealso \code{\link[grDevices]{col2rgb}} for translating \R colors to RGB vectors.
#' @export
#' @keywords internal
col2hexa <- function(x, alpha.channel=FALSE, alpha=1, charstring=TRUE){
  col <- col2rgb(x)
  if(alpha.channel){
    alpha <- as.integer(alpha*255)
    col <- rgb(red=col[1,], green=col[2,], blue=col[3,], alpha=alpha, maxColorValue = 255)
  } else {
    col <- rgb(red=col[1,], green=col[2,], blue=col[3,], maxColorValue = 255)
  }
  if(charstring){
    col <- paste("\"", col, "\"", sep="")
  }
  return(col)
}


#' Basemap Tiles Servers
#' 
#' Print a list of tiles servers ready-to-use with \code{\link{basemap}}.
#' 
#' @param print.servers logical. Should the names of the servers be printed?
#' @return Returns invisibly a matrix with servers names, urls and credits.
#' @export
bmSource <- function(print.servers=TRUE){
  bm.source <- matrix(c(
                        "mapquest.map",             "http://otile1.mqcdn.com/tiles/1.0.0/osm/{z}/{x}/{y}.jpg",
                        "mapquest.sat",             "http://otile1.mqcdn.com/tiles/1.0.0/sat/{z}/{x}/{y}.jpg",
                        "stamen.toner",             "http://stamen-tiles-{s}.a.ssl.fastly.net/toner/{z}/{x}/{y}.png",
                        "stamen.toner.hybrid",      "http://tile.stamen.com/toner-hybrid/{z}/{x}/{y}.png",
                        "stamen.toner.labels",      "http://tile.stamen.com/toner-labels/{z}/{x}/{y}.png",
                        "stamen.toner.lines",       "http://tile.stamen.com/toner-lines/{z}/{x}/{y}.png",
                        "stamen.toner.background",  "http://stamen-tiles-{s}.a.ssl.fastly.net/toner-background/{z}/{x}/{y}.png",
                        "stamen.toner.lite",        "http://stamen-tiles-{s}.a.ssl.fastly.net/toner-lite/{z}/{x}/{y}.png",
                        "stamen.watercolor",        "http://stamen-tiles-{s}.a.ssl.fastly.net/watercolor/{z}/{x}/{y}.png",
                        "cartodb.positron",         "http://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png",
                        "cartodb.positron.nolab",   "http://{s}.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png",
                        "cartodb.darkmatter",       "http://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png",
                        "cartodb.darkmatter.nolab", "http://{s}.basemaps.cartocdn.com/dark_nolabels/{z}/{x}/{y}.png"
                      ), ncol=2, byrow=T)
  
  mapquest.tiles.cr <- "Tiles: <a href=\"http://www.mapquest.com/\" target=\"_blank\" title=\"Tiles Courtesy of MapQuest\">MapQuest</a>"
  stamen.tiles.cr <- "Tiles: <a href=\"http://stamen.com\" title=\"Map tiles by Stamen Design, under CC BY 3.0.\">Stamen Design</a>"
  cartodb.tiles.cr <- "Tiles: <a href=\"http://cartodb.com/\" title=\"CartoDB\">CartoDB</a>"
  osm.data.cr <- "Data: <a href=\"http://openstreetmap.org\" title=\"Data by OpenStreetMap, under CC BY SA\">OSM</a>"
  nasa.data.cr <- "Data: <a href=\"\" title=\"NASA/JPL-Caltech and U.S. Depart. of Agriculture, Farm Service Agency\">NASA...</a>"
  
  
  mapquest.map.cr <- paste(mapquest.tiles.cr, osm.data.cr, sep=" | ")
  mapquest.sat.cr <- paste(mapquest.tiles.cr, nasa.data.cr, sep=" | ")
  stamen.cr <- paste(stamen.tiles.cr, osm.data.cr, sep=" | ")
  cartodb.cr <- paste(cartodb.tiles.cr, osm.data.cr, sep=" | ")
  vec.cr <- c(mapquest.map.cr, mapquest.sat.cr, rep(stamen.cr, 7), rep(cartodb.cr, 4))
  bm.source <- cbind(bm.source, vec.cr)
  if(print.servers){
    cat(paste(bm.source[,1], collapse="\n"))
  }
  invisible(bm.source)
}

#' Tiles Servers URL
#' 
#' Take a tiles server name (as returned by \code{\link{bmSource}}) and return its url.
#' 
#' @param x a character string of the name of the server.
#' @return The url of the server.
bmServer <- function(x){
  bm.source <- bmSource(print.servers=FALSE)
  res <- bm.source[bm.source[,1]==x,2]
  if(length(res) == 0L){
    res <- x
  }
  return(res)
}


#' Tiles Servers Attribution
#' 
#' Take a tiles server url and return its attribution.
#' 
#' @param x a character string of the url of the server.
#' @return The attribution of the server.
bmCredit <- function(x){
  bm.source <- bmSource(print.servers=FALSE)
  res <- bm.source[bm.source[,2]==x,3]
  if(length(res) == 0L){
    res <- NULL
  }
  return(res)
}


# Get the global bounding box of the data included in a map
.getExtBox <- function(sppts, spico, splns, sppol, spgrid){
  if(length(sppts) > 0){
    sppts.bb <- t(sapply(sppts, function(x) c(min(x$coords[, 1]), max(x$coords[, 1]),
                                     min(x$coords[, 2]), max(x$coords[, 2]))))
  } else {
    sppts.bb <- rep(NA, 4)
  }
  
  if(length(spico) > 0){
    spico.bb <- t(sapply(spico, function(x) c(min(x$coords[, 1]), max(x$coords[, 1]),
                                   min(x$coords[, 2]), max(x$coords[, 2]))))
  } else {
    spico.bb <- rep(NA, 4)
  }
  
  if(length(splns) > 0){
    splns.bb <- t(sapply(splns, function(x){z <- do.call("rbind", unlist(x$coords, recursive=FALSE));
                                     return(c(min(z[, 1]), max(z[, 1]), min(z[, 2]), max(z[, 2])))}))
  } else {
    splns.bb <- rep(NA, 4)
  }
  
  if(length(sppol) > 0){
    sppol.bb <- t(sapply(sppol, function(x){z <- do.call("rbind", unlist(x$coords, recursive=FALSE));
                              return(c(min(z[, 1]), max(z[, 1]), min(z[, 2]), max(z[, 2])))}))
  } else {
    sppol.bb <- rep(NA, 4)
  }
  
  if(length(spgrid) > 0){
    spgrid.bb <- t(sapply(spgrid, function(x) as.vector(x$x.bbox)[c(1,3,2,4)]))
  } else {
    spgrid.bb <- rep(NA, 4)
  }
  
  extbox <- rbind(sppts.bb, spico.bb, splns.bb, sppol.bb, spgrid.bb)
  extbox <- c(min(extbox[, 1], na.rm = TRUE),
              max(extbox[, 2], na.rm = TRUE),
              min(extbox[, 3], na.rm = TRUE),
              max(extbox[, 4], na.rm = TRUE))
  return(extbox)
}