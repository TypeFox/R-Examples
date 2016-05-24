#' Convert latlon data frame to gpx file.
#' 
#' Writes a gpx file of a dataframe with a call to \code{gpsbabel} in the shell
#' 
#' 
#' @param data Data frame with positions in columns \code{lat} and \code{lon}.
#' @param filename Name of gpx-file, defaults to 'tmp.gpx'.
#' @param type Type of gpx-file, one of \code{wpt} for waypoints, \code{rte}'
#' for route or \code{trk} for track.
#' @note Requires \code{gpsbabel} installation working from the command line.
#' @seealso \url{gpsbabel.org}
#' @keywords manip
#' @examples
#' 
#' \dontrun{
#' # some positions
#' pos <- rPeri(323)
#' frame2gpx(pos)
#' system("more tmp.gpx")
#' system("rm tmp.gpx")
#' }
#' 
#' @export frame2gpx
frame2gpx <-
function(data, filename = "tmp.gpx", type = "rte") {
  if (!is.data.frame(data)) stop("'data' not a dataframe")
  container <- tempfile("gpx")
  on.exit(unlink(container))
  colid <- match(c("lat", "lon"), names(data))
  data <- data[ , colid]
  data <- data.frame(data, 0:(nrow(data)-1))
  write.table(data, file = container,
    row.names = FALSE, col.names = FALSE, sep = ",")
  switch(type,
    wpt = system(paste("gpsbabel -i csv -f", container, 
      "-x transform -o gpx -F", filename)),
    rte = system(paste("gpsbabel -i csv -f", container, 
      "-x transform,rte=wpt,del -o gpx -F", filename)),
    trk = system(paste("gpsbabel -i csv -f", container, 
      "-x transform,trk=wpt,del -o gpx -F", filename))
  )
}

