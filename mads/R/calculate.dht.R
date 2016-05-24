#' Calculates the abundance for each species code including the unidentified 
#' codes if supplied.
#'  
#' @param species.name character vector of species codes
#' @param species.field.name character vector giving the field name of the 
#'   ddf data that contains the species codes
#' @param model.index named character vector which acts as a look up table for
#'   duplicate detection function models
#' @param ddf.results a list of ddf objects 
#' @param region.table dataframe of region records - Region.Label and Area
#' @param sample.table dataframe of sample records - Region.Label,
#'   Sample.Label, Effort
#' @param obs.table dataframe of observation records with fields object,
#'   Region.Label, and Sample.Label which give links to sample.table,
#'   region.table and the data records used in \code{model}
#' @param dht.options a list of the options to be supplied to mrds::dht
#' @return a list of dht objects, one for each species code
#' @author Laura Marshall
#' @seealso \code{mrds::dht}
#'
calculate.dht <- function(species.name, species.field.name, model.index, ddf.results, region.table, sample.table, obs.table, dht.options){
# calculate.dht function to calculate the abundance for each species code 
#
# Arguments:
#  species.name     -  character vector of species codes
#  model.index      -  named character vector which acts as a model look up table
#  ddf.results      -  a list of ddf objects
#  region.table     -  dataframe of region records
#  sample.table     -  dataframe of sample records
#  obs.table        -  dataframe of observation records
#
# Value:
#   list of dht objects
# 
# Function Calls: create.obs.table, mrds::dht 
#                                      
  #run dht
  dht.results <- list()
  for(sp in seq(along = species.name)){
    #get model name
    model.species.name <- model.index[[species.name[sp]]]
    #create obs.table [subset functionality in mrds relies on covariates in the ddf.data... this would involve changing stuff in Distance and/or mrds]
    if(!is.null(ddf.results[[model.species.name]])){
      obs.table.subset <- create.obs.table(obs.table, ddf.data = ddf.results[[model.species.name]]$data, subset.variable = species.field.name, subset.value = species.name[sp])
      #run dht analysis
      dht.results[[species.name[sp]]] <- try(dht(ddf.results[[model.species.name]], region.table, sample.table, obs.table.subset, options = dht.options))
      if(any(class(dht.results[[species.name[sp]]]) == "try-error")){
        warning("Error running dht for species ",species.name[sp], sep = "", fill = TRUE)
        dht.results[[species.name[sp]]] <- NULL
      }
    }else{
      dht.results[[species.name[sp]]] <- NULL
    }
  }
  return(dht.results) 
}

