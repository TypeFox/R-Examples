#' graticule: graticule lines for maps
#'
#' @docType package
#' @name graticule
NULL
limfun <- function(x, lim, nd = 60, meridian = TRUE) {
  ind <- 1:2
  if (!meridian) ind <- 2:1
  cbind(x = x, y = seq(lim[1], lim[2], length = nd))[, ind]
}

buildlines <- function(x) {
  do.call("rbind", lapply(seq_along(x), function(xx) {
    res <- data.frame(x[[xx]], rep(xx, nrow(x[[xx]])))
    names(res) <- c("x", "y", "id")
    res
  }))
}

## from raster findMethods("isLonLat")[["character"]]
# isLonLat <- function (x)
# {
#   res1 <- grep("longlat", as.character(x), fixed = TRUE)
#   res2 <- grep("lonlat", as.character(x), fixed = TRUE)
#   if (length(res1) == 0L && length(res2) == 0L) {
#     return(FALSE)
#   }
#   else {
#     return(TRUE)
#   }
# }

lonlatp4 <- function() {
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
}

#' Create graticule lines.
#'
#' Specify the creation of lines along meridians by specifying their placement
#' at particular \code{lons} (longitudes) and \code{lats} (latitudes) and their extents
#' with \code{xlim} (extent of parallel line in longitude) and \code{ylim} (extent of meridional line in latitude).
#'
#' Provide a valid PROJ.4 string to return the graticule lines in this projection. If this is not specified the graticule
#' lines are returned in their original longlat / WGS84.
#'
#' The arguments \code{xlim}, \code{ylim} and \code{nverts} are ignored if \code{tiles} is \code{TRUE}.
#' @param lons longitudes for meridional lines
#' @param lats latitudes for parallel lines
#' @param nverts number of discrete vertices for each segment
#' @param xlim maximum range of parallel lines
#' @param ylim maximum range of meridional lines
#' @param proj optional proj.4 string for output object
#' @param tiles if \code{TRUE} return polygons as output
#' @export
#'
#' @examples
#'\donttest{\dontrun{
#' library(rgdal)
#' x <- as.matrix(expand.grid(x = seq(100, 240, by = 15), y = seq(-85, -30, by = 15)))
#' prj <- "+proj=laea +lon_0=180 +lat_0=-70 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs "
#' px <- project(x, prj)
#' g <- graticule(unique(x[,1]), unique(x[,2]))
#' pg <- spTransform(g, CRS(prj))
#' plot(px, type = "n")
#' plot(pg, add = TRUE)
#'
#' g2 <- graticule(unique(x[,1]), unique(x[,2]), ylim = c(-90, 0), xlim = c(110, 250))
#' pg2 <- spTransform(g2, CRS(prj))
#' plot(px, type = "n")
#' plot(pg2, add = TRUE)
#'
#' prj <- "+proj=laea +lon_0=0 +lat_0=-90 +ellps=WGS84"
#' xx <- c(-120, -100, -80, -60, -40); yy <- c(-65, -55, -45)
#' g3 <- graticule(xx, yy, ylim = c(-70, -30), proj = prj)
#' g3labs <- graticule_labels(xx, c(-65, -45), xline = -85, yline = -30, proj = prj)
#' plot(g3)
#' text(g3labs, lab = parse(text = g3labs$lab))
#'
#' ## polygonal graticule on Orthographic projection
#' xx <- seq(-90, 90, length = 10) + 147
#' yy <- seq(-90, 90, length = 5)
#' portho <- "+proj=ortho +lon_0=147 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
#'  g <- graticule(xx, yy, proj = portho, tiles = TRUE)
#'  plot(g, col = c("black", "grey"))
#'
#'  library(maptools)
#'  data(wrld_simpl)
#'  w <- spTransform(subset(wrld_simpl, NAME == "Australia"), CRS(projection(g)))
#'  plot(w, add = TRUE, border = "dodgerblue")
#'  }}
#' @importFrom raster isLonLat raster rasterToPolygons extent values<- ncell
#' @importFrom sp SpatialLinesDataFrame Line Lines SpatialLines CRS spTransform
graticule <- function(lons, lats, nverts = 60, xlim, ylim, proj = NULL, tiles = FALSE) {
  if (is.null(proj)) proj <- lonlatp4()
  proj <- as.character(proj)  ## in case we are given CRS
  trans <- FALSE
  if (!raster::isLonLat(proj)) trans <- TRUE
  if (missing(lons)) {
    #usr <- par("usr")
    #if (all(usr == c(0, 1, 0, 1))) {
      lons <- seq(-180, 180, by = 15)
  }
  if (missing(lats)) {
      lats <- seq(-90, 90, by = 10)
  }
if (tiles) {
  ## build a raster and return the polygons
  rr <- raster::raster(extent(range(lons), range(lats)), nrows = length(lats)-1, ncols = length(lons)-1, crs = lonlatp4())
  values(rr) <- seq(ncell(rr))
  ## we need to not use raster for this
  nodes <- c(4, 8, 16)[pmax(1, findInterval(nverts, c(2, 3, 5)))]
  pp <- raster::rasterToPolygons(rr, dissolve = TRUE, n = nodes)
  if (trans) pp <- sp::spTransform(pp, sp::CRS(proj))
 return(pp)
}
  if (missing(xlim)) xlim <- range(lons)
  if (missing(ylim)) ylim <- range(lats)
  xline <- lapply(lons, limfun, lim = ylim, meridian = TRUE)
  yline <- lapply(lats, limfun, lim = xlim, meridian = FALSE)
  xs <- buildlines(xline)
  ys <- buildlines(yline)
  ys$id <- ys$id + max(xs$id)
  xs$type <- "meridian"
  ys$type <- "parallel"
  d <- rbind(xs, ys)
#  coordinates(d) <- ~x+y
 # proj4string(d) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  d0 <- split(d, d$id)
  l <- vector("list", length(d0))
  for (i in seq_along(d0)) l[[i]] <- sp::Lines(list(sp::Line(as.matrix(d0[[i]][, c("x", "y")]))), ID = as.character(i))
  l <- sp::SpatialLinesDataFrame(sp::SpatialLines(l, proj4string = sp::CRS(lonlatp4())),
                        data.frame(ID = as.character(seq_along(l))))
  if (trans) l <- sp::spTransform(l, sp::CRS(proj))
  l
  #d0$type <- c(rep("meridian", length(unique(xs$id))), rep("parallel", length(unique(ys$id))))

}

#' Create graticule labels.
#'
#' Returns a set of points with labels, for plotting in conjuction with \code{\link{graticule}}.
#'
#' SpatialPoints are returned in the projection of \code{proj} if given, or longlat / WGS84.
#' @param lons longitudes for meridional labels
#' @param lats latitudes for parallel labels
#' @param xline meridian/s for placement of parallel labels
#' @param yline parallel/s for placement of meridian labels
#' @param proj optional proj.4 string for output object
#' @export
#' @importFrom sp degreeLabelsEW degreeLabelsNS coordinates<- proj4string<-
#' @examples
#' xx <- c(100, 120, 160, 180)
#' yy <- c(-80,-70,-60, -50,-45, -30)
#' prj <- "+proj=lcc +lon_0=150 +lat_0=-80 +lat_1=-85 +lat_2=-75 +ellps=WGS84"
#' plot(graticule(lons = xx, lats = yy,  proj = prj))
#' labs <- graticule_labels(lons = xx, lats = yy, xline = 100, yline = -80,  proj = prj)
#' op <- par(xpd = NA)
#' text(labs, lab = parse(text = labs$lab), pos = c(2, 1)[labs$islon + 1], adj = 1.2)
#' par(op)
graticule_labels <- function(lons, lats, xline, yline, proj = NULL) {
  if (is.null(proj)) proj <- lonlatp4()
  proj <- as.character(proj)  ## in case we are given CRS
  trans <- FALSE
  if (!raster::isLonLat(proj)) trans <- TRUE
  if (missing(lons)) {
    #usr <- par("usr")
    #if (all(usr == c(0, 1, 0, 1))) {
    lons <- seq(-180, 180, by = 15)
  }
  if (missing(lats)) {
    lats <- seq(-90, 90, by = 10)
  }
  if (missing(xline)) xline <- lons
  if (missing(yline)) yline <- lats

  lonlabs <- expand.grid(x = lons, y = yline)
  lonlabs$lab <-  degreeLabelsEW(lonlabs$x)
  lonlabs$islon <- TRUE
  latlabs <- expand.grid(x = xline, y = lats)
  latlabs$lab <- degreeLabelsNS(latlabs$y)
  latlabs$islon <- FALSE
  l <- rbind(lonlabs, latlabs)
  coordinates(l) <- 1:2
  proj4string(l) <- CRS(lonlatp4())
  if (trans) {
    l <- sp::spTransform(l, CRS(proj))
  }
  l
}
