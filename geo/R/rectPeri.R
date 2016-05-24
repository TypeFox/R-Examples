#' Given rectangle code return perimeter as a polygon in lat lon
#' 
#' 
#' The outline/boundary of a statistical rectangle is returned as 5 positions,
#' the first and last of which are the same.
#' 
#' @name rectPeri
#' @aliases rPeri srPeri mrPeri drPeri
#' @param r,sr,mr,dr Rectangle codes.
#' @param dlat,dlon Dimensions of latitude and longitude given in minutes and
#' degrees for \code{mrPeri} and \code{drPeri}, respectively.
#' @return Rectangle outline as 5 positions.
#' @note Should perhaps be extended to give a list or dataframe of polygons for
#' more than one \code{r, sr, mr} or \code{dr}.
#' @seealso \code{\link{deg2rect}}, \code{\link{rectArea}},
#' \code{\link{geoarea}}.
#' @keywords arith manip
#' @examples
#' 
#'   geoplot(island, type = "n", grid = FALSE)
#'   geolines(rPeri(468))
#'   geolines(srPeri(4681))
#'

#' @export rPeri
#' @rdname rectPeri
rPeri <-
function(r)
{
  lat <- r2d(r)$lat
  lon <- r2d(r)$lon
  lat <- lat + c(1/4, 1/4,  - 1/4,  - 1/4, 1/4)
  lon <- lon + c(0.5,  - 0.5,  - 0.5, 0.5, 0.5)
  data.frame(lat = lat, lon = lon)
}

#' @export srPeri
#' @rdname rectPeri
srPeri <-
function(sr)
{
  lat <- sr2d(sr)$lat
  lon <- sr2d(sr)$lon
  lat <- lat + c(1/8, 1/8,  - 1/8,  - 1/8, 1/8)
  lon <- lon + c(0.25,  - 0.25,  - 0.25, 0.25, 0.25)
  data.frame(lat = lat, lon = lon)
}

#' @export mrPeri
#' @rdname rectPeri
mrPeri <-
function(mr, dlat = 5, dlon = 10)
{
  lat <- mr2d(mr, dlat = dlat, dlon = dlon)$lat
  lon <- mr2d(mr, dlat = dlat, dlon = dlon)$lon
  lat <- lat + c(dlat/120, dlat/120,  - dlat/120,  - dlat/120, dlat/120)
  lon <- lon + c(dlon/120,  - dlon/120,  - dlon/120, dlon/120, dlon/120)
  data.frame(lat = lat, lon = lon)
}

#' @export drPeri
#' @rdname rectPeri
drPeri <-
function(dr, dlat = 1, dlon = 2)
{
  lat <- dr2d(dr, dlat = dlat, dlon = dlon)$lat
  lon <- dr2d(dr, dlat = dlat, dlon = dlon)$lon
  lat <- lat + c(dlat/2, dlat/2,  - dlat/2,  - dlat/2, dlat/2)
  lon <- lon + c(dlon/2,  - dlon/2,  - dlon/2, dlon/2, dlon/2)
  data.frame(lat = lat, lon = lon)
}

