#' Given position return rectangle code.
#' 
#' 
#' Functions that convert positions in decimal degrees latitude and longitude
#' to rectangle codes for statistical rectangles, their subrectangles or
#' rectangle codes in systems of various resolutions as described below and
#' relating to (see \code{\link{rect2deg}}).
#' 
#' 
#' \itemize{
#' 
#' \item \code{r2d} with a resolution of 30 min latitue and 1 deg longitude
#' (the Icelandic numbering system, 'tilkynningaskyldureitir').
#' 
#' \item \code{sr2d} with a resolution of 15 min latitude by 30 min longitude
#' in the Icelandic numbering system for statistical rectangles which starts
#' counting at 60 deg N latitude, with sub-rectangles of 30 min lat by 1 deg
#' lon coded 1, 2, 3 and 4 for the NW, NA, SW and SA quadrants respectively
#' 
#' }
#' 
#' A small number (1e-06) is added to latitude and subtracted from longitude to
#' ensure rectangle membership of positions on border are \dQuote{logical} on
#' the nw-hempisphere.
#' 
#' @name deg2rect
#' @aliases d2r d2sr d2mr d2dr
#' @param lat,lon Position(s) as decimal degrees latitude and longitude.  If
#' \code{lat} is \code{list} its components \code{lat$lat} and \code{lat$lon}
#' are used for \code{lat} and \code{lon}.
#' @param dlat,dlon Rectangle height and width in degrees and minutes latitude
#' and longitude for \code{d2dr} and \code{d2mr} respectively.
#' @param startLat Starting latitude used in coding the rectangles.
#' @return Vector of rectangle codes in the chosen coding system.
#' @note These functions could be made hemisphere-aware, and no attention has
#' been given to making the functions work in the southern hemisphere, possibly
#' with the option of specifying starting latitudes.
#' @author
#' 
#' HB (\code{d2r, d2sr}, STJ (\code{d2mr, d2dr}).
#' @seealso
#' 
#' \code{\link{rect2deg}}
#' @keywords manip arith
#' @examples
#' 
#' 
#' ## tally positions in rectangles in object \code{\link{island}} giving
#' ## Iceland's coastline
#' 
#' data(island)
#' rects <- d2r(island)
#' table(rects)
#' 
#'

#' @export d2r
#' @rdname deg2rect
d2r <-
function(lat, lon = NULL)
{
	if(is.null(lon)) {
		lon <- lat$lon
		lat <- lat$lat
	}
	lat <- lat + 1e-06
	lon <- lon - 1e-06
	lon <-  - lon
	r <- (floor(lat) - 60) * 100 + floor(lon)
	ifelse(lat - floor(lat) > 0.5, r + 50, r)
}

#' @export d2sr
#' @rdname deg2rect
d2sr <-
function(lat, lon = NULL)
{
  if(is.null(lon)) {
    lon <- lat$lon
    lat <- lat$lat
  }
  lat <- lat + 1e-06
  lon <- lon - 1e-06
  lon <-  - lon
  r <- (floor(lat) - 60) * 100 + floor(lon)
  r <- ifelse(lat - floor(lat) > 0.5, r + 50, r)
  deg <- r2d(r)
  lon <-  - lon
  dlat <-  - (lat - deg$lat)
  dlon <-  - (lon - deg$lon)
  dl <- sign(dlat + 1e-07) + 2 * sign(dlon + 
    1e-07) + 4
  sr <- c(2, 0, 4, 0, 1, 0, 3)
  sr <- sr[dl]
  floor(r * 10 + sr)
}

#' @export d2mr
#' @rdname deg2rect
d2mr <-
function(lat, lon = NULL, dlat = 5, dlon = 10)
{
  if(is.null(lon)) {
    lon <- lat$lon
    lat <- lat$lat
  }
  lat <- lat + 1e-06
  lon <- lon - 1e-06
  lat <- geoconvert(lat, inverse = TRUE)
  lon <-  - geoconvert(lon, inverse = TRUE)
  mlat <- lat %% 10000 %/% 100
  mlon <- lon %% 10000 %/% 100
  mlat <- mlat %/% dlat
  mlon <- mlon %/% dlon
  lat <- lat %/% 10000
  lon <- lon %/% 10000
  lat * 1000000 + mlat * 10000 + lon * 100 + mlon
}

#' @export d2dr
#' @rdname deg2rect
d2dr <-
function(lat, lon = NULL, dlat = 1, dlon = 2, startLat = 50)
{
  if(is.null(lon)) {
    lon <- lat$lon
    lat <- lat$lat
  }
  lat <- lat + 1e-06
  lon <- lon - 1e-06
  hemi <- sign(lon)
  lat <- floor(lat)%%startLat
  lon <- floor(lon)
  hemi*(100*lat%/%dlat + hemi*floor(lon/dlon))
}

