#' Get all observations for a particular species or set of species
#'
#' @export
#' @param speciesid species id numbers, from 1 to infinity, potentially,
#'     use e.g., c(52, 53, etc.) if more than one species desired (numeric)
#' @param startdate start date of data period desired, see format in examples (character)
#' @param enddate end date of data period desired, see format in examples (character)
#' @template curl
#' @return A list with slots for taxa, stations, phenophase (metadata) and data
#' @examples \dontrun{
#' # Lookup names
#' lookup_names(name='Pinus', type='genus')
#'
#' # Get data on one species
#' npn_allobssp(speciesid = 52, startdate='2008-01-01', enddate='2011-12-31')
#'
#' # Get data on two species
#' npn_allobssp(speciesid = c(52, 53), startdate='2008-01-01', enddate='2011-12-31')
#'
#' # Get data on one species, convert to a single data.frame
#' npn_allobssp(speciesid = 52, startdate='2008-01-01', enddate='2011-12-31')
#' }

npn_allobssp <- function(speciesid, startdate = NULL, enddate = NULL, ...) {

  taxa <- taxonlist[as.numeric(as.character(taxonlist[,"species_id"])) %in% speciesid,c("species_id","genus","epithet","genus_epithet")]
  taxa$species_id <- as.numeric(as.character(taxa$species_id))

  args <- npnc(list(start_date = startdate, end_date = enddate))
  for (i in seq_along(speciesid)) {
    args[paste0('species_id[',i,']')] <- speciesid[i]
  }

  tt <- npn_GET(paste0(base(), 'observations/getAllObservationsForSpecies.json'), args, ...)

  #station_list <- data.frame(rbindlist(lapply(tt$station_list, data.frame), fill = TRUE, use.names = TRUE))
  phenophase_list <- lapply(tt$phenophase_list, function(x){ x[sapply(x, is.null)] <- "none"; x})
  phenophase_list <- data.frame(rbindlist(lapply(phenophase_list, data.frame), fill = TRUE, use.names = TRUE))

  statdf <- setDF(rbindlist(lapply(tt$station_list, function(z) {
    data.frame(t(unlist(pop(z, "species"))), stringsAsFactors = FALSE)
  }), fill = TRUE, use.names = TRUE))

  spp <- setNames(lapply(tt$station_list, function(z) {
    lapply(z$species, function(w) {
      df <- setDF(rbindlist(lapply(w, data.frame, stringsAsFactors = FALSE), fill = TRUE, use.names = TRUE))
      data.frame(id = names(w), df, stringsAsFactors = FALSE)
    })
  }), vapply(tt$station_list, "[[", 1, "station_id"))

  list(taxa = taxa, stations = statdf, phenophase = phenophase_list, data = spp)
}
