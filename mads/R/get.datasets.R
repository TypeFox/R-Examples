#' Extracts the data and models from the ddf objects
#'
#' Compiles a list of unique model combinations, which are stored as character 
#' vectors within the list. The associated data for the models in each of the
#' elements of this list are also obtained (these will only be one per set
#' of models as each model group must use the same data). It also constructs a 
#' look up table in the from of a named character vector to relate species code
#' to models as these have now been reduced to unique model combinations only.
#'  
#' @param model.names a list of character vectors of model names 
#'   with the elements named by species code
#' @param ddf.models a list of the ddf models named in model.names
#' @return list with the following elements:
#'   unique.model.names - a list of unique model combinations
#'   ddf.dat.master - a list of dataframes containing the data used to fit
#'     the unique models combinations defined in unique.model.names
#'   model.index - named character vector indicating which model is to be used
#'     for each species.                                                                       
#' @note Internal function not intended to be called by user.
#' @author Laura Marshall
#' @keywords data preparation
#'
get.datasets <- function(model.names, ddf.models){
# get.datasets function to extract unique models and corresponding data and provide the
# required look up table for those models used for more than one species.
#
# Arguments:
#   species.name a character vector of species codes
#   model.names a list of character vectors of model names for all species codes
#
# Value: list containing the following elements
#   unique.model.names a list of unique character vectors of model names 
#   ddf.dat.master a list of dataframes containing the data used to fit the models named in unique.model.names
#   model.index look up table in the form of a named character vector
#
# Function Calls: none
#
  species.name <- names(model.names)
  model.id <- species.name[1]
  unique.ddf.models <- list()
  #Add the models for the first species as they must be unique
  unique.ddf.models[[species.name[1]]] <- sort(model.names[[species.name[1]]]) 
  #For each species code (except the first)
  for(sp in seq(along = species.name)[-1]){
    #find mode names for new species code
    current.model.names <- sort(model.names[[species.name[sp]]])  
    #set up a boolean flag to check if the models are already in the unique.ddf.models list 
    exist <- FALSE
    for(uniq in seq(along = unique.ddf.models)){
      #check to see if the combination of models are already in the unique.ddf.models list
      if(length(unique.ddf.models[[uniq]]) == length(current.model.names)){
        #if there are the same number of models compare
        compare <- unique.ddf.models[[uniq]] == current.model.names
      }else{       
        #otherwise they must be different
        compare <- FALSE
      }
      #If they are the same length check if they are the same
      if(length(which(!compare)) == 0){
        #If they are the same store these models in common.species
        exist <- TRUE
        common.species <- names(unique.ddf.models)[uniq]
      }
    }
    #If the models are already in the unique.ddf.models
    if(exist){   
      #Update model.id
      model.id[sp] <- common.species
    }else{
      #Otherwise add them to unique.ddf.models
      unique.ddf.models[[species.name[sp]]] <- sort(model.names[[species.name[sp]]]) 
      #Update model.id
      model.id[sp] <- species.name[sp]  
    }
  }  
  names(model.id) <- species.name
  #Extract data from models in unique.ddf.models
  ddf.dat.master <- list()
  unique.ddf.model.names <- names(unique.ddf.models)
  for(sp in seq(along = unique.ddf.model.names)){
    #Only take the first model as all models have the same data [checked above]
    ddf.dat.master[[unique.ddf.model.names[sp]]] <- ddf.models[[unique.ddf.models[[sp]][1]]]$data  
  } 
  return(list(ddf.dat.master = ddf.dat.master, unique.model.names = unique.ddf.models, model.index = model.id))
}