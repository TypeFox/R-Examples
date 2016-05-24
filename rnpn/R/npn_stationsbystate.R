#' Get number of stations by state.
#'
#' @export
#' @template curl
#' @return Number of stations by state as a data.frame.
#' @examples \dontrun{
#' head( npn_stationsbystate() )
#' }
npn_stationsbystate <- function(...) {
  tt <- npn_GET(paste0(base(), 'stations/getStationCountByState.json'), list(), ...)
  states <- sapply(tt, function(x){
    if (is.null(x[[1]]) == TRUE) {
      x[[1]] <- "emptyvalue"
    } else{
      x[[1]] <- x[[1]]
    }
  })
  data <- sapply(tt, "[[", "number_stations")
  structure(data.frame(states, data, stringsAsFactors = FALSE), .Names = c("state", "number_stations"))
}
