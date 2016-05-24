#google static API plotting

gmap.bbox <- function(lon, lat, zoom, wdthpx, htpx) {
  pixels <- 2^zoom*256
  centerX <- (lon+180)/360.0*pixels
  centerY <- (180-sm.lat2y(lat))/360.0*pixels

  swX <- centerX-wdthpx/2
  swlon <- swX/pixels*360-180
  neX <- centerX+wdthpx/2
  nelon <- neX/pixels*360-180

  swY <- centerY-htpx/2
  swlat <- -sm.y2lat(swY/pixels*360-180)
  neY <- centerY+htpx/2
  nelat <- -sm.y2lat(neY/pixels*360-180)

  matrix(c(swlon, swlat, nelon, nelat), ncol=2, byrow=FALSE)
}



google.getimage <- function(maptype, lon, lat, zoom, wdpx,
                            htpx, scale=1, key=NULL, cachedir=NULL, forcedownload=FALSE) {
  if(is.null(key)) {
    key <- "AIzaSyA7TenP6BraUbiUvea_bhJqdeaK1nyOSjE"
  }
  format="png"

  url <- paste0("https://maps.googleapis.com/maps/api/staticmap?maptype=", "satellite",
                "&center=", lat, ",", lon, "&zoom=", zoom,
                "&size=", wdpx, "x", htpx, "&key=", key, "&format=", format,
                "&scale=", scale)
  folder <- tile.cachedir("google", cachedir)
  filename <- paste0(digest::digest(url), ".", format)
  tofile <- file.path(folder, filename)

  if(!file.exists(tofile) || forcedownload) {
    message("Downloading to ", tofile)
    utils::download.file(url, tofile, quiet=TRUE)
  }

  png::readPNG(tofile)
}

gmap.plot <- function(bbox, maptype="satellite", forcedownload=FALSE,
                           cachedir=NULL, res=150, project=TRUE, key=NULL, ...) {

  if(project) {
    epsg <- 3857
  } else {
    epsg <- 4326
  }
  #setup plot
  bboxplot <- .projectbbox(bbox, epsg)
  coords <- sp::coordinates(t(bboxplot))
  spoints = sp::SpatialPoints(coords, proj4string = sp::CRS(paste0("+init=epsg:", epsg)))
  plotargs <- list(...)
  if(is.null(plotargs$xlim))
    xlim <- bboxplot[1,]
  if(is.null(plotargs$ylim))
    ylim <- bboxplot[2,]
  sp::plot(spoints, pch=".", xlim=xlim, ylim=ylim, ...)

  #get extents
  ext <- graphics::par("usr")
  fullareabbox <- matrix(c(ext[1], ext[3], ext[2], ext[4]), ncol=2, byrow=FALSE)
  fullareabboxll <- .revprojectbbox(fullareabbox, epsg)

  bboxgoog <- sm.projectbbox(fullareabboxll)
  if((bboxgoog[1,2] - bboxgoog[1,1]) < 0) {
    bboxgoog[1,2] <- bboxgoog[1,2]+360
  }
  midx <- (bboxgoog[1]+bboxgoog[3])/2
  midy <- (bboxgoog[2]+bboxgoog[4])/2
  midlatlon <- sm.tolatlon(midx, midy)

  zoom <- tile.autozoom(res, epsg)
  #y / x
  aspect <- (bboxgoog[2,2] - bboxgoog[2,1]) / (bboxgoog[1,2] - bboxgoog[1,1])

  anglewidth <- fullareabboxll[1,2] - fullareabboxll[1,1]
  if(anglewidth < 0) {
    anglewidth <- anglewidth+360
  }
  totpixels <- 2^zoom*256
  widthpx <- (anglewidth / 360.0) * totpixels
  heightpx <- widthpx * aspect

  #calc scale factor based on res parameter
  widthin <- graphics::grconvertX(ext[2], from="user", to="inches") -
    graphics::grconvertX(ext[1], from="user", to="inches")
  actualres <- widthpx / widthin
  if(actualres < res) {
    scale <- 2
  } else {
    scale <- 1
  }

  cat("aspect:", aspect, " widthpx:", widthpx, " heightpx:", heightpx, " zoom:", zoom,
      " center:", midlatlon[1], " ", midlatlon[2], "\n")
  cat("bbox query:", bboxgoog, "\n")
  bboximgll <- gmap.bbox(midlatlon[1], midlatlon[2], zoom, widthpx, heightpx)
  cat("bbox return: ", sm.projectbbox(bboximgll), "\nbbox reg:", bboximgll, "\n")

  image <- google.getimage(maptype, midlatlon[1], midlatlon[2], zoom, round(widthpx),
                           round(heightpx), scale=scale, key=key, cachedir=cachedir,
                           forcedownload=forcedownload)
  tile.plotarray(image, .projectbbox(bboximgll, epsg))
}
