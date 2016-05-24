#' Plot Squares in Different Colors
#' 
#' Plot Icelandic statistical squares in different colors, based on specified
#' levels.
#' 
#' If \code{levels} is a named vector, the names are used as level-specific
#' colors, ignoring the last element.
#' 
#' The \code{levels} should go lower and higher than the range of \code{z}
#' values.
#' 
#' @param sr squares or subsquares.
#' @param z values for each square.
#' @param levels threshold levels for splitting the values into categories.
#' @param grid passed to \code{geoplot}.
#' @param \dots passed to \code{geoplot}.
#' @return Invisible \code{NULL}.
#' @note Convenience wrapper for the \code{reitaplott} function.
#' 
#' The 3-digit Icelandic statistical squares are 0.5 degrees latitude and 1
#' degree longitude, defined in the range 60--70\eqn{^{\circ}}{°}N and
#' 0--50\eqn{^{\circ}}{°}W.
#' 
#' The bottom-left corner of square 0 is at 60\eqn{^{\circ}}{°}N
#' 0\eqn{^{\circ}}{°}. From there, the numbering system is best explained with
#' a table:
#' 
#' \tabular{lrrrrr}{ \tab 49\eqn{^{\circ}}{°}W \tab 48\eqn{^{\circ}}{°}W \tab
#' \dots{} \tab 1\eqn{^{\circ}}{°}W \tab 0\eqn{^{\circ}}{°}W \cr
#' 69.5\eqn{^{\circ}}{°}N \tab \samp{999} \tab \samp{998} \tab \dots{} \tab
#' \samp{951} \tab \samp{950} \cr 69.0\eqn{^{\circ}}{°}N \tab \samp{949} \tab
#' \samp{948} \tab \dots{} \tab \samp{901} \tab \samp{900} \cr \dots{} \tab
#' \dots{} \tab \dots{} \tab \dots{} \tab \dots{} \tab \dots{} \cr
#' 60.5\eqn{^{\circ}}{°}N \tab \samp{099} \tab \samp{098} \tab \dots{} \tab
#' \samp{051} \tab \samp{050} \cr 60.0\eqn{^{\circ}}{°}N \tab \samp{049} \tab
#' \samp{048} \tab \dots{} \tab \samp{001} \tab \samp{000} \cr }
#' 
#' The 4-digit subsquares divide each square into quadrants (\acronym{NW, NE,
#' SW, SE}), appending a number to indicate the quadrant:
#' 
#' \tabular{ll}{ \samp{***1} \tab \samp{***2}\cr \samp{***3} \tab
#' \samp{***4}\cr }
#' @author Arni Magnusson.
#' @seealso \code{\link{geoplot}} and \code{\link{reitaplott}} are
#' the underlying drawing functions.
#' 
#' \code{\link{colorRampPalette}} can be used to create a ramp of
#' level-specific colors.
#' @keywords hplot spatial
#' @examples
#' 
#' # Use existing palette
#' geoSR(561:560, c(0.3,3), levels=c(0,1,10))
#' 
#' # Pass colors along with levels
#' lev <- c(0, 1, 10)
#' names(lev) <- c("brown", "orange", NA)
#' geoSR(561:560, c(0.3,3), lev)
#' 
#' # Subsquares
#' geoSR(5611:5612, c(0.3,3), lev)
#' 
#' # Color ramp
#' z <- (0:10) / 10
#' lev <- (0:10) / 10
#' ramp <- colorRampPalette(c("khaki1","gold","orange","darkorange2","red",
#'                            "darkred","black"))
#' names(lev) <- ramp(length(lev))
#' geoSR(724:715, z, lev)
#' 
#' @export geoSR
geoSR <- function(sr, z, levels, grid=FALSE, ...)
{
  if(!is.null(names(levels)))
  {
    opal <- palette(c("black", names(levels)))
    on.exit(palette(opal))
  }

  geoplot(grid=grid, ...)
  invisible(capture.output(reitaplott(reitur=sr, smareitur=NULL, z=z, levels=levels, colors=seq(levels), density=0)))
  geopolygon(geo::island)
  geolines(geo::island)
}
