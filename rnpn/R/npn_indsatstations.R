#' Get all observations for a particular species or set of species.
#'
#' @export
#' @param stationid Required. Use e.g., c(4881, 4882, etc.) if more than one species desired
#' (numeric)
#' @template curl
#' @return Observations for each species by date in a data.frame
#' @examples \dontrun{
#' npn_indsatstations(stationid = c(507, 523))
#' }
npn_indsatstations <- function(stationid, ...) {
  args <- list()
  for (i in seq_along(stationid)) {
    args[paste0('station_ids[',i,']')] <- stationid[i]
  }
  tt <- npn_GET(paste0(base(), 'individuals/getIndividualsAtStations.json'), args, ...)
  ldfply(tt)
}
