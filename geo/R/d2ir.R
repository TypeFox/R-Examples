#' Convert between Geographic Coordinates and ICES Rectangles
#' 
#' This function converts geographic coordinates to ICES North Atlantic
#' statistical rectangles and rectangle codes to rectangle center coordinates.
#' 
#' The default \code{useI=FALSE} is in accordance with the prescription in ICES
#' CM77/Gen:3, but \code{useI=TRUE} has been done on occasion.
#' 
#' @aliases d2ir ir2d
#' @param lat vector of latitudes, or a list containing \code{lat} and
#' \code{lon}.
#' @param lon vector of longitudes (ignored if \code{lat} is a list).
#' @param useI whether to use the letter \sQuote{\samp{I}} in statistical
#' rectangle codes.
#' @param ir ICES rectangle code, e.g. \samp{37F3}
#' @return Vector of strings containing ICES statistical rectangle codes or
#' conversely, center coordinates of rectangles which codes were given as
#' input.
#' @note The statistical rectangles are coded as follows:
#' 
#' The latitudinal rows are 30' high, numbered from \samp{01} at the southern
#' boundary of the ICES statistical area (36°00'N, see ICES No. 2 chart)
#' northwards to \samp{99}. The northern limit of the statistical recangle
#' system is thus 85°30'N. The longitudinal columns are 1° wide, coded in an
#' alphanumeric system which starts at the western boundary of the ICES
#' statistical area (44°00'W, see ICES No.1 chart) with \samp{A0}, continuing
#' \samp{A1}, \samp{A2}, \samp{A3} to 40°W. East of 40°W the system continues
#' \samp{B0}, \samp{B1}, \samp{B2}, \ldots{}, \samp{B9}, \samp{C0}, \samp{C1},
#' \samp{C2}, \ldots{}, \samp{M8}, using a different letter for each 10° block
#' and covering the entire west-east extent of the ICES statistical area, thus:
#' 
#' \tabular{lrcrl}{ \samp{A0} \tab 44°00'W \tab -- \tab 40°00'W \tab
#' \samp{A3}\cr \samp{B0} \tab 40°00'W \tab -- \tab 30°00'W \tab \samp{B9}\cr
#' \samp{C0} \tab 30°00'W \tab -- \tab 20°00'W \tab \samp{C9}\cr \samp{D0} \tab
#' 20°00'W \tab -- \tab 10°00'W \tab \samp{D9}\cr \samp{E0} \tab 10°00'W \tab
#' -- \tab 0°00' \tab \samp{E9}\cr \samp{F0} \tab 0°00' \tab -- \tab 10°00'E
#' \tab \samp{F9}\cr \samp{G0} \tab 10°00'E \tab -- \tab 20°00'E \tab
#' \samp{G9}\cr \samp{H0} \tab 20°00'E \tab -- \tab 30°00'E \tab \samp{H9}\cr
#' \samp{J0} \tab 30°00'E \tab -- \tab 40°00'E \tab \samp{J9}\cr \samp{K0} \tab
#' 40°00'E \tab -- \tab 50°00'E \tab \samp{K9}\cr \samp{L0} \tab 50°00'E \tab
#' -- \tab 60°00'E \tab \samp{L9}\cr \samp{M0} \tab 60°00'E \tab -- \tab
#' 68°30'E \tab \samp{M8}\cr }
#' 
#' Note that letter \sQuote{\samp{I}} is omitted, unless \code{useI=TRUE}.
#' 
#' When designating a statistical rectangle, the north coordinate is stated
#' first. Thus, the recangle of which the southwest corner is 54°00'N, 3°00'E
#' is coded \samp{37F3}.
#' @author Sigurdur Thor Jonsson with contributions 
#' from Asta Gudmundsdottir and Arni Magnusson.
#' @seealso \code{\link{d2r}} and \code{\link{r2d}} convert between geographic
#' coordinates and Icelandic rectangles (a local coding system using the same
#' rectangle size as the ICES coding system).
#' @references ICES C.M. 1977/Gen:3.
#' @keywords manip arith
#' @examples
#' 
#' d2ir(54.25, 3.5)
#' d2ir(c(50,60), c(-20,-10))
#' 
#' ir2d(d2ir(54, 3))
#' ## center positions for bottom left and approx top right rects
#' ir2d("01A0")
#' ir2d("99M7")
#' ir2d(c("01A0","99M7"))
#' ## note that ICES CM1977/Gen:3 indicates half-size rects on eastern margin!
#'
#' @export d2ir
d2ir <-
function(lat, lon = NULL, useI = FALSE)
{
  if(is.null(lon)) {
    lon <- lat$lon
    lat <- lat$lat
  }
  lat <- lat + 1e-06
  lon <- lon + 1e-06
  outside <- lat < 36 | lat >= 85.5 | lon <= -44 | lon > 68.5
  if(any(outside)) warning("Positions outside of ICES statistical area")
  lat <- floor(lat * 2) - 71
  lat <- ifelse(lat < 10, paste("0", lat, sep = ""), lat)
  if(useI) lettersUsed <- LETTERS[1:12]
  else lettersUsed <- LETTERS[c(1:8,10:13)]
  lon1 <- lettersUsed[(lon + 60) %/% 10]
  lon2 <- ifelse(lon1 == "A", floor(lon %% 4), floor(lon %% 10))
  ir <- paste(lat, lon1, lon2, sep = "")
  ir[outside] <- NA
  ir
}
#' @export ir2d
#' @rdname d2ir
ir2d <-
function(ir, useI = FALSE)
{
  lat <- substring(ir, 1, 2)
  lat <- as.numeric(lat)
  lat <- (lat +71)/2 + 0.25
  lon1 <- substring(ir, 3, 3)
  lon1 <- toupper(lon1)
  lon1 <- match(lon1, LETTERS)
  if(!useI) lon1 <- ifelse(lon1 > 8, lon1 - 1, lon1)
  lon1 <- lon1 - 2
  lon2 <- substring(ir, 4)
  lon2 <- as.numeric(lon2)
  lon <- ifelse(lon1 < 0,
                -44 + lon2 + 0.5,
                -40 + 10*lon1 + lon2 + 0.5)
  data.frame(lat = lat, lon = lon)
}
