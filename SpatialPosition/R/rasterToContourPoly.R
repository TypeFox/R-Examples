#' @title Create a SpatialPolygonsDataFrame from a Raster
#' @name rasterToContourPoly
#' @description 
#' This function creates a contour SpatialPolygonsDataFrame from a raster.
#' @param r raster; the raster must contain only positive values.
#' @param nclass numeric; a number of class.
#' @param breaks numeric; a vector of break values. 
#' @param mask SpatialPolygonsDataFrame; mask used to clip contour shapes. 
#' The mask should have a smaller extent than r.
#' @return The ouput of the function is a SpatialPolygonsDataFrame. 
#' The data frame of the outputed SpatialPolygonsDataFrame contains four fields: 
#' id (id of each polygon), min and max (minimum and maximum breaks of the polygon), 
#' center (central values of classes)
#' @details This function uses the rgeos package.
#' @seealso \link{stewart}, \link{rasterStewart}, \link{plotStewart}, 
#' \link{quickStewart}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @import sp
#' @import raster
#' @examples
#' data("spatData")
#' \dontrun{
#' mystewart <- stewart(knownpts = spatPts, varname = "Capacite",
#'                      typefct = "exponential", span = 1000, beta = 3,
#'                      resolution = 50, longlat = FALSE,
#'                      mask = spatMask)
#' # Create a raster of potentials values
#' mystewartraster <- rasterStewart(x = mystewart)
#' # Create contour SpatialLinesDataFrame
#' contourpoly <- rasterToContourPoly(r = mystewartraster,
#'                                    nclass = 6,
#'                                    mask = spatMask)
#' # Created breaks
#' bks <- unique(c(contourpoly$min, contourpoly$max))
#' # Display the map
#' library(cartography)
#' opar <- par(mar = c(0,0,1.2,0))
#' choroLayer(spdf = contourpoly,
#'            df = contourpoly@data,
#'            var = "center", legend.pos = "topleft",
#'            breaks = bks, border = "grey90",
#'            lwd = 0.2,
#'            legend.title.txt = "Potential number\nof beds in the\nneighbourhood",
#'            legend.values.rnd = 0)
#' plot(spatMask, add = TRUE)
#' propSymbolsLayer(spdf = spatPts, df = spatPts@data, var = "Capacite",
#'                  legend.title.txt = "Number of beds",
#'                  col = "#ff000020")
#' layoutLayer(title = "Global Accessibility to Public Hospitals",
#'             south = TRUE, sources = "", author = "")
#' par(opar)
#' }
#' @export
rasterToContourPoly <- function(r, nclass = 8, breaks = NULL, mask = NULL){
  if (!requireNamespace("rgeos", quietly = TRUE)) {
    stop("'rgeos' package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if(!'package:rgeos' %in% search()){
    attachNamespace('rgeos')
  }
  
  # get initial min and max values
  rmin <- raster::cellStats(r, min, na.rm = TRUE)
  rmax <- raster::cellStats(r, max, na.rm = TRUE)
  
  # default breaks and nclass
  if(is.null(breaks)){
    breaks <- seq(from = rmin, 
                  to = rmax, 
                  length.out = (nclass+1))
  }else{
    breaks <- c(rmin, breaks[breaks > rmin & breaks < rmax], rmax)
    breaks <- unique(breaks)
    breaks <- sort(breaks)
  }
  
  # get raster resolution and projection  
  myres <- raster::res(r)[1]
  myproj <- sp::CRS(sp::proj4string(r))
  
  # extent or mask the raster around the mask
  if (is.null(mask)){
    mask <- masker(r)
    maskbuff <- rgeos::gBuffer(mask, byid = FALSE, width = 5 * myres )
    r <- raster::extend(r, maskbuff, value=-1)
  }else{
    mask <- rgeos::gUnaryUnion(mask)
    maskbuff <- rgeos::gBuffer(mask, byid = FALSE, width = 5 * myres )
    r <- raster::mask(r, maskbuff, updatevalue = -1)  
    # test mask extent
    if(rgeos::gWithin(masker(r), mask)){stop("mask should be smaller than r",
                                             call. = FALSE)}
  }
  
  # adjust breaks if necessary
  rmin <- min(r[r!=-1])
  rmax <- max(r[r!=-1])
  breaks <- c(rmin, breaks[breaks > rmin & breaks < rmax], rmax)
  breaks <- unique(breaks)
  breaks <- sort(breaks)
  finalBreaks <- breaks
  # zero level problem
  if(breaks[1] <= 0){
    zv <- TRUE
    breaks <- breaks + 1
    r <- r + 1
  }else{
    zv <- FALSE
  }
  nclass <- length(breaks)-1
  breaks <- breaks[-(nclass+1)]
  
  # deal with NAs
  r[is.na(r)] <- 0
  
  # test breaks
  if(length(breaks)<2){stop("breaks values do not fit the raster values", 
                            call. = FALSE)}
  
  # build the contour lines polygones
  cl <- rasterToContour(r, levels = breaks)
  cl$level <- as.numeric(as.character(cl$level))
  SPlist <- list()
  SPlevels <- character()
  for (i in cl$level){ 
    linex <- cl[cl@data$level == i,]
    linex <- linex@lines
    linex <- linex[[1]]
    linex <- linex@Lines
    Plist <- NULL
    Plist <- list()
    for (j in 1:length(linex)){
      x <- linex[[j]]@coords
      x <- sp::Polygon(coords =  x, hole = F)
      x <- sp::Polygons(srl = list(x), ID = j)
      Plist[j] <- x
    }  
    x <- sp::SpatialPolygons(Srl = Plist)
    x <- rgeos::union(x = x)
    if (class(x) != "SpatialPolygonsDataFrame"){
      x <- sp::SpatialPolygonsDataFrame(Sr = x, 
                                        data = data.frame(
                                          level = rep(i, length(x))))
    } else {
      x <- x[x@data$count < 2,]
      x@data <- data.frame(level = rep(i, dim(x)[1]))
    }
    SPlist <- c(SPlist , x@polygons  )
    SPlevels <- c(SPlevels,x@data$level)
  }
  for (i in 1:length(SPlist)){
    SPlist[[i]]@ID <- as.character(i)
  }
  x <- sp::SpatialPolygons(Srl = SPlist, proj4string = myproj)
  x <- sp::SpatialPolygonsDataFrame(Sr = x, 
                                    data = data.frame(levels = SPlevels))
  
  bks <- data.frame(b =c(breaks, rmax), t = finalBreaks)
  
  # manage attributes data of the contour spdf
  x@data <- data.frame(id = paste("id_",row.names(x),sep=""),
                       min = bks[match(x$levels, bks[,1]),2], 
                       max = bks[match(x$levels, bks[,1])+1,2],
                       center = NA, 
                       stringsAsFactors = FALSE)
  x$center <- (x$min+x$max) / 2 
  row.names(x) <- x$id
  
  # clip the contour spdf with the mask
  final <- rgeos::gIntersection(spgeom1 = x, spgeom2 = mask, byid = TRUE, 
                                id = row.names(x))
  df <- data.frame(id = sapply(methods::slot(final, "polygons"), 
                               methods::slot, "ID"))
  row.names(df) <- df$id
  final <- sp::SpatialPolygonsDataFrame(Sr = final, data = df)
  final@data <- data.frame(id = final$id, x[match(final$id, x$id),2:4])
  final@plotOrder <- 1:nrow(final)
  return(final)
}


masker <- function(r){
  xy <- sp::coordinates(r)[which(!is.na(values(r))),]
  i <- grDevices::chull(xy)
  b <- xy[c(i,i[1]),]
  mask <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(b, hole = FALSE)), 
                                                ID = "1")), 
                              proj4string = sp::CRS(sp::proj4string(r)))
  return(mask)
}
