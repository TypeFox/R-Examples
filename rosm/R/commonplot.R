#functions used by both google and osm

tile.cachedir <- function(type, cachedir=NULL) {
  if(is.null(cachedir)) {
    cachedir <- "rosm.cache"
  }
  folder <- file.path(cachedir, type)
  created <- dir.create(folder, showWarnings=FALSE, recursive=TRUE)
  folder
}

tile.plotarray <- function(image, box) {
  graphics::rasterImage(image, box[1,1], box[2,1], box[1,2], box[2,2])
}

tile.autozoom <- function(res=150, epsg=4326) {
  ext <- graphics::par("usr")
  midy <- mean(c(ext[3], ext[4]))
  rightmid <- .tolatlon(ext[2], midy, epsg)
  centermid <- .tolatlon(mean(c(ext[1], ext[2])), midy, epsg)
  leftmid <- .tolatlon(ext[1], midy, epsg)

  anglewidth1 <- rightmid[1] - centermid[1]
  if(anglewidth1 < 0) {
    anglewidth1 <- anglewidth1+360
  }

  anglewidth2 <- rightmid[1] - centermid[1]
  if(anglewidth2 < 0) {
    anglewidth2 <- anglewidth2+360
  }
  anglewidth <- anglewidth1+anglewidth2

  #PROBLEMS WITH WIDE EXTENTS LIKE THE WORLD
  widthin <- graphics::grconvertX(ext[2], from="user", to="inches") -
    graphics::grconvertX(ext[1], from="user", to="inches")
  widthpx <- widthin * res

  zoom = log2((360.0 / anglewidth) * (widthpx / 256.0))

  as.integer(floor(zoom))
}
