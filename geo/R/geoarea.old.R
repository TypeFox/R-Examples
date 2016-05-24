#' Old method for calculating geographical area
#' 
#' Calculates the area inside a geographical region by splitting it into a
#' grid, determining which gridpoints are inside the region and then finds the
#' area.
#' 
#' 
#' @param reg Region of interest, list with components \code{lat} and
#' \code{lon}.
#' @param n Numeber of gridpoints to use for calculating the area.
#' @param robust Robust or not, defaults to \code{TRUE}.
#' @return Geographical area in square nautical miles.
#' @note Difference between old and new \code{geoarea}-method might be
#' explained better, if indeed it is necessary to keep the old method.
#' @seealso \code{\link{geoarea}.}
#' @keywords arith
#' @export geoarea.old
"geoarea.old"<-
function(reg, n,robust=T) {
  reg$lat <- (reg$lat * pi)/180
  reg$lon <- (reg$lon * pi)/180
  rlat <- range(reg$lat[!is.na(reg$lat)])
  dlat <- (rlat[2] - rlat[1])
  rlon <- range(reg$lon[!is.na(reg$lon)])
  dlon <- (rlon[2] - rlon[1]) * cos((rlat[2] + rlat[1])/2)
  ratio <- (dlat/dlon)
  nlat <- floor(sqrt(n) * sqrt(ratio))
  nlon <- round(sqrt(n)/sqrt(ratio))
  dlon <- dlon/cos((rlat[2] + rlat[1])/2)
  lat <- rlat[1] + (c(0:(nlat - 1)) * dlat)/nlat + dlat/(2 * nlat)
  lon <- rlon[1] + (c(0:(nlon - 1)) * dlon)/nlon + dlon/(2 * nlon)
  darea <- (lon[2] - lon[1]) * (lat[2] - lat[1]) * 40528473
  latgr <- c(matrix(lat, nlat, nlon))
  longr <- c(t(matrix(lon, nlon, nlat)))
  area <- dlon * dlat * 40528473 * cos((rlat[2] + rlat[1])/2)
  border <- adapt(reg$lat, reg$lon)
  inside<- rep(0, length(latgr))
  if(robust) {
    a <- a1 <- rep(0, length(reg$lat))
    
    inside<- .C("marghc",
	       as.double(longr),
	       as.double(latgr),
	       as.integer(length(latgr)),
	       as.double(border$lon),
	       as.double(border$lat),
	       as.integer(length(border$lat)),
	       as.integer(border$lxv),
	       as.integer(length(border$lxv)),
	       as.integer(inside),
	       as.double(a),
	       as.double(a1))
    inside<- inside[[9]]
  }
  else {
    tmpinni <- rep(0,length(border$lxv))
    inside<- .C("geomarghc",
	       as.double(longr),
	       as.double(latgr),
	       as.integer(length(latgr)),
	       as.double(border$lon),
	       as.double(border$lat),
	       as.integer(border$lxv),
	       as.integer(length(border$lxv)),
	       as.integer(inside),
	       as.integer(tmpinni))
    inside<- inside[[8]]
  }
  ind <- c(1:length(inside))
  ind <- ind[inside!= 0]
  mlat <- mean(latgr[ind])
  mlon <- mean(longr[ind])
  cmlat <- mean(cos(latgr[ind]))
  cl <- mean(cos(latgr))
  rat <- length(ind)/length(inside)	# fraction outside
  inside.area <- (rat * area * cmlat)/cl
  return(inside.area)
}

