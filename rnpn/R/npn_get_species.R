#' Get scientific names.
#'
#' @export
#' @param ids One or more ITIS taxonomic serial numbers (tsn), or NPN ID numbers.
#' @param state Required. A US state, two-letter abbreviation.
#' @param kingdom Optional. A taxonomic kingdom.
#' @param genus A genus name
#' @param species A specific epithet, the second part of a full species name
#' @param name A common name
#' @param network The primary key of the network for which to filter species.
#' @param year Year of obseration
#' @param groups One or more primary keys associated with a species type.
#' @param stationid Station ID. Use e.g., c(4881, 4882, etc.) if more than one species desired
#' @template curl
#' @return data.frame of species and their IDs
#' @examples \dontrun{
#' head( npn_species() )
#' npn_species_itis(ids = 27806)
#' npn_species_itis(ids = c(27806,36616))
#' npn_species_id(ids = 3)
#' npn_species_state(state = "HI")
#' npn_species_state(state = "HI", kingdom = "Plantae")
#' npn_species_sci(genus = "Clintonia", species = "borealis")
#' npn_species_comm(name = "thickleaved wild strawberry")
#' npn_species_comm(name = c("thickleaved wild strawberry","bluebead"))
#' npn_species_search(groups = 3, year = 2010)
#' npn_species_search(groups = c(3,9), year = 2010)
#'
#' library('httr')
#' npn_species_itis(ids = 27806, config=verbose())
#' }
npn_species <- function(...) {
  ldfply(npn_GET(paste0(base(), 'species/getSpecies.json'), list(), ...))
}

#' @export
#' @rdname npn_species
npn_species_itis <- function(ids, ...) {
  tt <- lapply(ids, function(z){
    npn_GET(paste0(base(), 'species/getSpeciesByItis.json'), list(itis_sn = z), ...)
  })
  ldfply(tt)
}

#' @export
#' @rdname npn_species
npn_species_id <- function(ids, ...) {
  tt <- lapply(ids, function(z){
    npn_GET(paste0(base(), 'species/getSpeciesById.json'), list(species_id = z), ...)
  })
  ldfply(tt)
}

#' @export
#' @rdname npn_species
npn_species_state <- function(state, kingdom = NULL, ...) {
  args <- npnc(list(state = state, kingdom = kingdom))
  ldfply(npn_GET(paste0(base(), 'species/getSpeciesByState.json'), args, ...))
}

#' @export
#' @rdname npn_species
npn_species_sci <- function(genus, species, ...) {
  args <- list(genus = genus, species = species)
  data.frame(npn_GET(paste0(base(), 'species/getSpeciesByScientificName.json'), args, ...),
             stringsAsFactors = FALSE)
}

#' @export
#' @rdname npn_species
npn_species_comm <- function(name, ...) {
  tt <- lapply(name, function(z){
    npn_GET(paste0(base(), 'species/getSpeciesByCommonName.json'), list(common_name = z), ...)
  })
  ldfply(tt)
}

#' @export
#' @rdname npn_species
npn_species_search <- function(network=NULL, year=NULL, groups=NULL, stationid=NULL, ...) {
  args <- npnc(list(network_id = network, observation_year = year))
  for (i in seq_along(groups)) {
    args[paste0('group_ids[',i,']')] <- groups[i]
  }
  for (i in seq_along(stationid)) {
    args[paste0('station_ids[',i,']')] <- stationid[i]
  }

  ldfply(npn_GET(paste0(base(), 'species/getSpeciesFilter.json'), args, ...))
}
