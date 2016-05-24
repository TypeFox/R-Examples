#'@title get lake names in GLTC dataset
#'@description Lake names that are part of the Global Lake Temperature Collaboration can be 
#'accessed and returned. See associated publication and references therein for 
#'details including units and data provenance. 
#'@return a character vector of valid lake names from 
#'the Global Lake Temperature Collaboration dataset.
#'@seealso \code{\link{get_climate_names}}, \code{\link{get_metadata_names}}
#'@examples
#'get_lake_names()
#'@export
get_lake_names <- function(){

  lake_names <- unique(gltc_metadata$Lake.name)
  
  return(lake_names)
}