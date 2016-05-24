#' avisSpeciesId
#' 
#' Returns the id of the selected species
#' 
#' @usage avisSpeciesId(nameraw)
#' @param nameraw scientific name of the species (e.g. "Pica pica")
#' @return an integer
#' @export 
#' @examples \dontrun{
#' avisSpeciesId("Pica pica")
#' }
#' 

avisSpeciesId <- function (nameraw)
{
  if(!avisHasSpecies(nameraw)){
    stop(paste("Species '", nameraw,"' not found in proyectoavis.com"))
  }

  allspecies <- avisAllSpecies()

  return (as.integer(allspecies[.avisNormalizeSpeciesName(nameraw)]))
}


#' avisHasSpecies
#' 
#' check if a species name exists in Proyecto AVIS. 
#' 
#' @usage avisHasSpecies(nameraw)
#' @param nameraw scientific name of the species (e.g. "Pica pica")
#' @return Logical: returns TRUE for species with observations in the database and 
#' FALSE otherwise
#' @examples \dontrun{
#' avisHasSpecies("Pica pica")
#' avisHasSpecies("Pica pic")
#' }
#' @export
#' 
avisHasSpecies <- function (nameraw)
{
  name <- .avisNormalizeSpeciesName(nameraw)
  allspecies <- avisAllSpecies()

  return (is.element(name, names(allspecies)))
}

#' avisAllSpecies
#' 
#' Returns a vector with the ids of the species in Proyecto AVIS
#' 
#' @usage avisAllSpecies()
#' @note This functions does not allow arguments
#' @return returns a vector
#' @export 
#' @examples \dontrun{
#' avisAllSpecies()
#' }
#' 
#' 
avisAllSpecies <- function()
{
  .avisCacheReturnOrSetup(".ravis_species_id_list", ".avisGetServerEspecies")
}

# fetches species list from server
.avisGetServerEspecies <- function()
{
  .avisVerboseMessage("INFO: fetching species list from proyectoavis.com server")

  species<- .avisApiBusOrden()

  return (species)
}


#' avisSpeciesSummary
#' 
#' Download a table with a summary of the records 
#' stored in Proyecto AVIS (http://proyectoavis.com) 
#' aggregated by species; number of observations of 
#' each species, number of individuals recorded, 
#' number of different UTMs (10x10km) with observations, 
#' number of birdwatchers that recorded the species
#' 
#' @usage avisSpeciesSummary()
#' @note This functions does not allow arguments
#' @return returns a dataframe
#' @export 
#' @examples \dontrun{
#' avis_summary<- avisSpeciesSummary()
#' # general overview of the data aggregated by species
#' par (mfrow =c(2,2))
#' hist (avis_summary$Observations, col="red", border=FALSE, main=NULL)
#' hist (avis_summary$Individuals, col="red", border=FALSE, main=NULL)
#' hist (avis_summary$UTM.10x10, col="red", border=FALSE, main=NULL)
#' hist (avis_summary$Birdwatchers, col="red", border=FALSE, main=NULL)
#' }
#' 

avisSpeciesSummary <- function ()
{
  .avisCacheReturnOrSetup(".ravis_species_summary", function(){

    # fetches form server
    .avisVerboseMessage("INFO: fetching species summary from proyectoavis.com server")

    spsummary<- .avisApiBusEspecie()

    return(spsummary)
  });
}

.avisNormalizeSpeciesName<- function(raw)
{
  return (tolower(raw))
}
