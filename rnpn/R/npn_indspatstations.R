#' Get all observations for a particular species or set of species.
#'
#' @export
#' @param speciesid Required. Species id numbers, from 1 to infinity, potentially,
#'     use e.g., c(52, 53, etc.) if more than one species desired (numeric)
#' @param stationid Required. Use e.g., c(4881, 4882, etc.) if more than one species desired
#' (numeric)
#' @param year Year (numeric).
#' @template curl
#' @return Observations for each species by date.
#' @examples \dontrun{
#' npn_indspatstations(speciesid = 35, stationid = c(60, 259), year = 2009)
#' npn_indspatstations(35, c(60, 259), 2009)
#' }

npn_indspatstations <-  function(speciesid, stationid, year = NULL, ...) {
  args <- npnc(list(year = year))
  for (i in seq_along(speciesid)) {
    args[paste0('species_id[',i,']')] <- speciesid[i]
  }
  for (i in seq_along(stationid)) {
    args[paste0('station_ids[',i,']')] <- stationid[i]
  }
  tt <- npn_GET(paste0(base(), 'individuals/getIndividualsOfSpeciesAtStations.json'), args, ...)
  ldfply(tt)
}
