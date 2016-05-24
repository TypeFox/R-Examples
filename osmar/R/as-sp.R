#' @include as-osmar.R
{}



#' Convert osmar object to sp object
#'
#' Convert an osmar object to a \link[sp]{sp} object.
#'
#' @param obj An \code{\link{osmar}} object
#' @param what A string describing the sp-object; see Details section
#' @param crs A valid \code{\link[sp]{CRS}} object; default value is
#'   given by \code{\link{osm_crs}}-function
#' @param simplify Should the result list be simplified to one element
#'   if possible?
#'
#' @details
#'   Depending on the strings given in \code{what} the
#'   \code{\link{osmar}} object will be converted into in a list of
#'   objects given by the \link[sp]{sp}-package:
#'
#'   \describe{
#'
#'     \item{\code{what = "points"}}{the object will be converted
#'       in a \code{\link[sp]{SpatialPointsDataFrame}}. The data slot is
#'       filled with the attrs slot of \code{obj$nodes}.}
#'
#'     \item{\code{what = "lines"}}{the object will be converted in
#'       a \code{\link[sp]{SpatialLinesDataFrame}}. It is build with all
#'       possible elements which are in \code{obj$ways}
#'       \code{obj$relations}. The data slot is filled with elements
#'       of both.}
#'
#'     \item{\code{what = "polygons"}}{the object will be converted
#'       in a \code{\link[sp]{SpatialPolygonsDataFrame}}. It consists of
#'       elements which are in \code{obj$ways} slot.}
#'
#'  }
#'
#'  Every conversion needs at least a non-empty
#'  \code{obj$nodes$attrs}-slot because spatial information are stored
#'  in there.
#'
#' @return
#'   A list of one or more \link[sp]{sp} objects; see Details section.
#'
#' @examples
#'   data("muc", package = "osmar")
#'   muc_points <- as_sp(muc, "points")
#'   muc_lines <- as_sp(muc, "lines")
#'   muc_polygons <- as_sp(muc, "polygons")
#'
#'   bbox(muc_points)
#'
#' @export
as_sp <- function(obj, what = c("points", "lines", "polygons"),
                  crs = osm_crs(), simplify = TRUE) {

  stopifnot(require("sp"))
  stopifnot(is_osmar(obj))
  stopifnot(any(has_data(obj)))

  what <- match.arg(what, several.ok = TRUE)

  ret <- lapply(what,
                function(w) {
                  fun <- sprintf("as_sp_%s", w)
                  do.call(fun, list(obj, crs))
                })
  names(ret) <- what

  if ( length(ret) == 1 & simplify )
    ret <- ret[[1]]

  ret
}




#' CRS for OpenStreetMap
#'
#' Coordinate Reference System used in OpenStreetMap.
#'
#' @param crs A valid proj4 string
#'
#' @details
#'   The default value is the WGS84 Ellipsoid which is used in GPS,
#'   therefore it is used in OpenStreetMap.
#'
#' @return
#'   A \code{\link[sp]{CRS}} object
#'
#' @examples
#'   osm_crs()
#'   class(osm_crs())
#'
#' @export
osm_crs <- function(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") {
  stopifnot(require("sp"))
  ret <- CRS(crs)
  ret
}



### Internal converter functions: ####################################


### ... to SpatialPoints:

as_sp_points <- function(obj, crs = osm_crs()) {
  ## TODO: possibility for multipoints-object (point-relations e.g.)

  fullcheck <- has_data(obj)
  if( !fullcheck["nodes"] ) {
    warning("no nodes")
    return(NULL)
  }

  dat <- unique(obj$nodes$attrs)
  coords <- data.frame(lon = dat$lon, lat = dat$lat, row.names = dat$id)
  ret <- SpatialPointsDataFrame(coords = coords, proj4string = crs,
                                data = dat, match.ID = "id")
  ret
}



### ... to SpatialLines:

as_sp_lines <- function(obj, crs = osm_crs()){
  fullcheck <- has_data(obj)
  if ( !fullcheck["nodes"] ) {
    warning("no nodes")
    return(NULL)
  }
  if ( !fullcheck["ways"] ) {
    warning("no ways")
    return(NULL)
  }

  way_ids <- unique(obj$ways$refs$id)
  way_lns <- vector("list", length(way_ids))
  for(i in 1:length(way_lns)) {
    way_lns[[i]]<- Lines(ways_nodes2Line(way_ids[i], obj$ways, obj$nodes),
                         way_ids[i])
  }

  if(fullcheck[3] == FALSE) {
    return(make_SLDF(obj, way_lns, crs, "way"))
  }

  rel_ids <- unique(obj$relations$refs$id)
  rel_lns <- vector("list", length(rel_ids))
  for(i in 1:length(rel_ids)) {
    rel_lns[[i]] <- Lines(rels_ways_nodes2Line(rel_ids[i], obj$relations,
                                               obj$ways, obj$nodes), rel_ids[i])
  }

  ret <- make_SLDF(obj, c(way_lns, rel_lns), crs, "relation")

  ret
}



rels_ways_nodes2Line <- function(relID, rels, ways, nodes){
  #ref <- subset(rels$refs, id == relID)  # CMD check note: no visible binding
  ref <- rels$refs[rels$refs$id == relID, ]
  #wayref <- subset(ref, type == "way")$ref  # CMD check note: no visible binding
  wayref <- ref[ref$type == "way", ]$ref
  wayln <-lapply(wayref, "ways_nodes2Line", ways, nodes)
#  relref<- subset(ref, type=="relation")$ref
#  falls ways der relations noch eingebaut werden sollen
  wayln <- wayln[!sapply(wayln, is.null)]
  wayln
}



ways_nodes2Line <- function(wayID, ways, nodes){
  #nds <- subset(ways$refs, id==wayID)$ref  # CMD check note: no visible binding
  nds <- ways$refs[ways$refs$id == wayID, ]$ref
  if ( length(nds) == 0) {
    return(NULL)
  }
  geo <- nodes$attrs[match(nds,nodes$attrs$id), c("lon","lat")]
  if(sum(is.na(geo)==TRUE)>=1)
    geo<- geo[!is.na(geo[,1]),]
  ret <- Line(geo)
  ret
}



make_SLDF <- function(obj, lns, crs, what){
  lns <- remove_emptyLines(lns)
  splns <- SpatialLines(lns, proj4string = crs)
  if ( what == "way") {
    dat <- cbind(obj$ways$attrs, type = as.factor("way"))
  }
  if ( what == "relation" ) {
    dat <- rbind(cbind(obj$ways$attrs,type = as.factor("way")),
                 cbind(obj$relations$attrs,type = as.factor("relation")))
  }
  ret <- SpatialLinesDataFrame(splns, data = dat, match.ID = "id")
  ret
}



remove_emptyLines <- function(LINES) {
  LINES[sapply(1:length(LINES), function(k) length(LINES[[k]]@Lines))!=0]
}



### ...  to SpatialPolygons:

as_sp_polygons <- function(obj, crs = osm_crs()){
  fullcheck <- has_data(obj)
  if ( !fullcheck["nodes"] ) {
    warning("no nodes")
    return(NULL)
  }
  if ( !fullcheck["ways"] ) {
    warning("no ways")
    return(NULL)
  }

  way_ids <- unique(obj$ways$refs$id)
  way_pols <- vector("list", length(way_ids))
  for(i in 1:length(way_pols)) {
    way_pols[[i]]<- ways_nodes2Polygon(way_ids[i], obj$ways, obj$nodes)
    if ( class(way_pols[[i]]) == "Polygon" )
      way_pols[[i]]<- Polygons(list(way_pols[[i]]), way_ids[i])
  }
  polys_position<- which(!sapply(way_pols, is.list))
  way_pols <- way_pols[polys_position]

  if( length(way_pols) == 0 ) {
    warning("no polygon-like objects in \"ways\"")
    return(NULL)
  }
  ## relations don't have many polygonlike objects

  spols <- SpatialPolygons(way_pols, proj4string=crs)
  dat <- obj$ways$attrs[polys_position,]
  ret <- SpatialPolygonsDataFrame(spols, data=dat, match.ID="id")
  ret
}



ways_nodes2Polygon <- function(wayID, ways, nodes){
  #nds <- subset(ways$refs, id==wayID)$ref  # CMD check note: no visible binding
  nds <- ways$refs[ways$refs$id == wayID, ]$ref
  if(length(nds)==0)
    return(list(NULL))

  geo <- nodes$attrs[match(nds,nodes$attrs$id), c("lon","lat")]
  if(sum(is.na(geo)==TRUE)>=1)
    geo<- geo[!is.na(geo[,1]),]
  if(sum(is_poly(geo)) != 2)
    return(list(NULL))

  ret <- Polygon(geo)
  ret
}



is_poly <- function(matrix){
  matrix[1,] == matrix[nrow(matrix), ]
}

