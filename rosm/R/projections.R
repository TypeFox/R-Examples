#projection functions

.tolatlon <- function(x, y, epsg=NULL, projection=NULL) {
  rgdal::CRSargs(sp::CRS("+init=epsg:3857")) #hack to load rgdal namespace
  if(is.null(epsg) && is.null(projection)) {
    stop("epsg and projection both null...nothing to project")
  } else if(!is.null(epsg) && !is.null(projection)) {
    stop("epsg and projection both specified...ambiguous call")
  }

  if(is.null(projection)) {
    projection <- sp::CRS(paste0("+init=epsg:", epsg))
  }

  coords <- sp::coordinates(matrix(c(x,y), byrow=TRUE, ncol=2))
  spoints <- sp::SpatialPoints(coords, projection)
  spnew <- sp::spTransform(spoints, sp::CRS("+init=epsg:4326"))
  c(sp::coordinates(spnew)[1], sp::coordinates(spnew)[2])
}

.fromlatlon <- function(lon, lat, epsg=NULL, projection=NULL) {
  rgdal::CRSargs(sp::CRS("+init=epsg:3857")) #hack to load rgdal namespace
  if(is.null(epsg) && is.null(projection)) {
    stop("epsg and projection both null...nothing to project")
  } else if(!is.null(epsg) && !is.null(projection)) {
    stop("epsg and projection both specified...ambiguous call")
  }

  if(is.null(projection)) {
    projection <- sp::CRS(paste0("+init=epsg:", epsg))
  }

  coords <- sp::coordinates(matrix(c(lon,lat), byrow=TRUE, ncol=2))
  spoints <- sp::SpatialPoints(coords, sp::CRS("+init=epsg:4326"))
  spnew <- sp::spTransform(spoints, projection)
  c(sp::coordinates(spnew)[1], sp::coordinates(spnew)[2])
}

.projectbbox <- function(bbox, toepsg=NULL, projection=NULL) {
  rgdal::CRSargs(sp::CRS("+init=epsg:3857")) #hack to load rgdal namespace
  if(is.null(toepsg) && is.null(projection)) {
    stop("toepsg and projection both null...nothing to project")
  } else if(!is.null(toepsg) && !is.null(projection)) {
    stop("toepsg and projection both specified...ambiguous call")
  }

  if(is.null(projection)) {
    projection <- sp::CRS(paste0("+init=epsg:", toepsg))
  }
  coords <- sp::coordinates(t(bbox))
  spoints = sp::SpatialPoints(coords, proj4string = sp::CRS("+init=epsg:4326"))
  newpoints <- sp::spTransform(spoints, projection)
  newbbox <- t(sp::coordinates(newpoints))

  if(newbbox[1,1] > newbbox[1,2]) { #if min>max
    maxx <- .fromlatlon(180, bbox[2, 1], projection=projection)[1]
    newbbox[1,1] <- newbbox[1,1]-maxx*2
  }
  newbbox
}

.revprojectbbox <- function(bbox, fromepsg=NULL, projection=NULL) {
  rgdal::CRSargs(sp::CRS("+init=epsg:3857")) #hack to load rgdal namespace
  if(is.null(fromepsg) && is.null(projection)) {
    stop("fromepsg and projection both null...nothing to project")
  } else if(!is.null(fromepsg) && !is.null(projection)) {
    stop("fromepsg and projection both specified...ambiguous call")
  }
  if(is.null(projection)) {
    projection <- sp::CRS(paste0("+init=epsg:", fromepsg))
  }
  coords <- sp::coordinates(t(bbox))
  spoints = sp::SpatialPoints(coords, proj4string = projection)
  newpoints <- sp::spTransform(spoints, sp::CRS("+init=epsg:4326"))
  t(sp::coordinates(newpoints))
}


sm.y2lat <- function(a) {
  return(180.0/pi * (2 * atan(exp(a*pi/180)) - pi/2))
}

sm.lat2y <- function(a) {
  return(180.0/pi * log(tan(pi/4+a*(pi/180)/2)))
}

sm.projectbbox <- function(bbox) {
  bbox[2,] <- sm.lat2y(bbox[2,])
  bbox
}

sm.revprojectbbox <- function(bbox) {
  bbox[2,] <- sm.y2lat(bbox[2,])
  bbox
}

sm.tolatlon <- function(x, y) {
  c(x, sm.y2lat(y))
}
