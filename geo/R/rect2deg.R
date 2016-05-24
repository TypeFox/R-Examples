#' Given rectangle code return its center position.
#' 
#' Functions that convert statistical rectangle codes under: 1) a traditional
#' Icelandic system ('tilkynningaskyldurreitakerfid') (see \code{\link{d2r}}
#' and \code{\link{d2sr}}) and 2) set up in systems based on minutes and
#' degrees (see \code{\link{d2mr}} and \code{\link{d2dr}}) to decimal
#' representation of rectangles center positions in degreees latitude and
#' longitude.
#' 
#' \itemize{
#' 
#' \item \code{r2d} with a resolution of 30 min latitue and 1 deg longitude
#' (the Icelandic numbering system, 'tilkynningaskyldureitir').
#' 
#' \item \code{sr2d} with a resolution of 15 min latitude by 30 min longitude
#' in the Icelandic numbering system for statistical rectangles which starts
#' counting at 60 deg N latitude, with sub-rectangles of 30 min lat by 1 deg
#' lon coded 1, 2, 3 and 4 for the NW, NA, SW and SA quadrants respectively,
#' 
#' \item \code{mr2d} with resolution given in \code{dlat} by \code{dlon}
#' minutes lat and lon,
#' 
#' \item \code{dr2d} for rectangles with resolution given in \code{dlat} by
#' \code{dlon} degrees lat and lon, code system starting at latitude
#' \code{startLat}.
#' 
#' }
#' 
#' @name rect2deg
#' @aliases r2d sr2d mr2d dr2d
#' @param r Rectangle code \code{r} in the 'tillkynningaskyldu-system', e.g
#' from \code{\link{deg2rect}}.
#' @param sr Rectangle code \code{sr} for subrectangle in
#' 'tilkynningaskyldu-system', e.g. from \code{\link{deg2rect}}.
#' @param mr Rectangle code \code{mr} based on minutes, e.g. from
#' \code{\link{deg2rect}}.
#' @param dr Rectangle code \code{dr} based on degrees, e.g. from
#' \code{\link{deg2rect}}.
#' @param dlat Rectangle height in minutes or degrees latitude for \code{mr2d}
#' and \code{dr2d} respectively.
#' @param dlon As \code{dlat} except now width in longitude.
#' @param startLat Starting latitude for coding the rectangles.
#' @return dataframe of center positions (latitude \code{lat} and longitude
#' \code{lon}) of rectangles in one of the coding systems
#' @note Mostly used for plotting.
#' @author
#' 
#' HB (for \code{r, sr}), STJ (for \code{mr, dr}).
#' @seealso \code{\link{deg2rect}}
#' @keywords arith manip
#' @examples
#' 
#'   r2d(d2r(lat = 65 + 1/4, lon = -19 - 1/2))
#'   d2r(r2d(519))
#'

#' @export r2d
#' @rdname rect2deg
r2d <-
function(r)
{
	lat <- floor(r/100)
	lon <- (r - lat * 100) %% 50
	halfb <- (r - 100 * lat - lon)/100
	lon <-  - (lon + 0.5)
	lat <- lat + 60 + halfb + 0.25
	data.frame(lat = lat, lon = lon)
}

#' @export sr2d
#' @rdname rect2deg
sr2d <-
function(sr)
{
  r <- floor(sr/10)
  sr <- sr - r * 10
  lat <- floor(r/100)
  lon <- (r - lat * 100) %% 50
  halfb <- (r - 100 * lat - lon)/100
  lon <-  - (lon + 0.5)
  lat <- lat + 60 + halfb + 0.25
  l1.lat <- c(0, 0.125, 0.125, -0.125, -0.125)
  l1.lon <- c(0, -0.25, 0.25, -0.25, 0.25)
  lat <- lat + l1.lat[sr + 1]
  lon <- lon + l1.lon[sr + 1]
  data.frame(lat = lat, lon = lon)
}

#' @export mr2d
#' @rdname rect2deg
mr2d <-
function(mr, dlat = 5, dlon = 10)
{
  lat <- mr %/% 1000000
  mr <- mr %% 1000000
  mlat <- mr %/% 10000
  mr <- mr %% 10000
  lon <- mr %/% 100
  mlon <- mr %% 100
  lat <- 10000 * lat + 100 * (dlat * mlat + dlat/2)
  lon <- 10000 * lon + 100 * (dlon * mlon + dlon/2)
  data.frame(lat = geoconvert(lat), lon =  - geoconvert(lon))
}

#' @export dr2d
#' @rdname rect2deg
dr2d <-
function(dr, dlat = 1, dlon = 2, startLat = 50)
{
  hemi <- sign(dr + 1e-06)
  lat <- startLat + dlat*(abs(dr)%/%100)
  lon <- dlon*(dr%%(hemi*100))
  data.frame(lat = lat + dlat/2, lon =  lon + dlon/2)
}

