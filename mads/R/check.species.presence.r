#' Checks the list of species presence definitions supplied by the user
#'
#' Performs a numer of checks to ensure the list has been defined correctly.
#' If no list was specified one is created assuming all species are present in 
#' all strata.
#' 
#' @param species.presence a list with an element for each  strata
#'   which contains a vector of corresponding species codes for those
#'   species which occur in that strata or NULL if it is to be data driven.
#' @param species.name a vector of species names for which model names were 
#'   supplied in the ddf.models list passed to execute.multi.analysis by the 
#'   user.
#' @param strata.name character vector of strata names
#' @return updated species.presence list 
#' @note Internal function not intended to be called by user.
#' @author Laura Marshall
#' @seealso \code{execute.multi.analysis}
#' @keywords input validation
#'
check.species.presence <- function(species.presence, species.name, strata.name){
# 
# check.species.presence function to check the list of species code definitions supplied by the user for each strata 
#
# Arguments:  #
#  species.cpresence - a list of vectors containing the species codes which correspond to the 
#    list element name
#  species.name             - a vector of species names obtained from the list of models
#
# Value:
#   the updated species.code.definitions list
#
# Function Calls: none 
# 
  #create if null
  if(is.null(species.presence)){ 
    species.precence <- list()
    for(st in seq(along = strata.name)){
      species.presence[[strata.name[st]]] <- species.name
    }  
  }
  return(species.presence)     
}



