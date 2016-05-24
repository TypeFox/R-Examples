#' Checks the list of species code deifinitions supplied by the user
#'
#' Firstly, it checks that there are not multiple definitions for the same species code.
#' Secondly, it checks that the names of the ddf models correspond to the names 
#' in the species code definitions. Lastly, it checks that definitions are 
#' provided for all species codes and if not it adds the missing definitions 
#' assuming they correspond to identified categories
#' 
#' @param species.code.definitions a list with an element for each 
#'   unidentified code which contains a vector of corresponding identified 
#'   species codes or NULL if not required
#' @param species.name a vector of species names for which model names were 
#'   supplied in the ddf.models list passed to execute.multi.analysis by the 
#'   user.
#' @return updated species.code.definitions list 
#' @note Internal function not intended to be called by user.
#' @author Laura Marshall
#' @seealso \code{execute.multi.analysis}
#' @keywords input validation
#'
check.species.code.definitions <- function(species.code.definitions, species.name){
# 
# check.species.code.definitions function to check the list of species code deifinitions supplied by the user 
#
# Arguments:  #
#  species.code.definitions - a list of vectors containing the species codes which correspond to the 
#    list element name
#  species.name             - a vector of species names obtained from the list of models
#
# Value:
#   the updated species.code.definitions list
#
# Function Calls: none 
# 
  #check to see if there are unidentified species codes
  unidentified = FALSE
  for(i in seq(along = species.code.definitions)){
    if(is.null(species.code.definitions[[i]])){      
      stop(paste("No species codes specified for ",names(species.code.definitions)[i]," in the species code definitions list.", sep = ""), call. = FALSE) 
    }else if(length(species.code.definitions[[i]]) > 1){
      unidentified = TRUE
      #check to see that unidentified codes are only being prorated to identified codes.
      for(j in seq(along = species.code.definitions[[i]])){
        if(species.code.definitions[[i]][j] == names(species.code.definitions)[i]){          
          stop("Incorrect species code definition for species ",names(species.code.definitions)[i],". Unidentified code cannot be prorated to itself.", call. = FALSE)
        }else if(species.code.definitions[[i]][j]%in%names(species.code.definitions)){   
          if(length(species.code.definitions[[species.code.definitions[[i]][j]]]) > 1){           
            stop("Incorrect species code definition for species ",names(species.code.definitions)[i],". An Unidentified code cannot be prorated to another unidentified code.", call. = FALSE)
          }
        }
      } 
    }else if(length(species.code.definitions[[i]]) == 1 & species.code.definitions[[i]] != names(species.code.definitions)[i]){   
      stop("Incorrect species code definition for species ",names(species.code.definitions)[i],". If only a single code is entered it must match the name of the list element.", call. = FALSE)      
    }else if(length(species.code.definitions[[i]]) == 0){ 
      
      stop(paste("No species codes specified for ",names(species.code.definitions)[i]," in the species code definitions list.", sep = ""), call. = FALSE) 
    }
  } 
  #get all unique codes in the list    
  all.codes <- definition.names <- names(species.code.definitions)
  for(sp in seq(along = species.code.definitions)){
    all.codes <- c(all.codes, species.code.definitions[[sp]])
  }
  unique.codes <- unique(all.codes)
  compare.codes <- ifelse(length(species.name) == length(unique.codes), sort(species.name) == sort(unique.codes), NA)
  #check that there are not multiple definitions for the same species code
  if(length(species.code.definitions) != length(unique(definition.names))){   
    stop("Multiple species code entries in the species code definitions list.", call. = FALSE)
  #check that the names of the ddf models correspond to the names in the species code definitions  
  }else if(is.na(compare.codes) | !is.na(compare.codes) & length(which(compare.codes)) == length(species.name)){   
    stop("Species mismatch in ddf models and species code definitions. Models not suppled for all species or models supplied for species not included in species code definitions.", call. = FALSE)
  #if there are no problems and definitions are provided for all species return the list unchanged
  }
  if(length(species.code.definitions) == length(unique.codes)){
    return(list(unidentified = unidentified, species.code.definitions = species.code.definitions)) 
  #if there are missing definitions add them in assuming they are identified categories
  }else {     
    for(sp in seq(along = species.name)){
      #check if there is an entry
      if(is.na(match(species.name[sp], definition.names))){
        #if not add one
        species.code.definitions[[species.name[sp]]] <- c(species.name[sp])
      }      
    }
    return(list(unidentified = unidentified, species.code.definitions = species.code.definitions))
  } 
}



