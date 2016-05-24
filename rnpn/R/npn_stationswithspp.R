#' Get a list of all stations which have an individual whom is a member of a
#'    set of species.
#'
#' @export
#' @param speciesid Required. Species id numbers, from 1 to infinity, potentially,
#'    use e.g., c(52, 53, etc.) if more than one species desired (numeric)
#' @template curl
#' @return Stations' latitude and longitude, names, and ids.
#' @examples \dontrun{
#' npn_stationswithspp(speciesid = c(52,53,54))
#' npn_stationswithspp(speciesid = 53)
#' }

npn_stationswithspp <- function(speciesid, ...) {
  args <- list()
  for (i in seq_along(speciesid)) {
    args[paste0('species_id[',i,']')] <- speciesid[i]
  }
  ldfply(npn_GET(paste0(base(), 'stations/getStationsWithSpecies.json'), args, ...))
}
