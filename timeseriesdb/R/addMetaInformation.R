#' Add Meta Information to R Environments 
#' 
#' This function adds meta information to environments that 
#' are explicitly meant to store Meta Information. This function 
#' can be used separately in interactive R Session or to facilitate
#' mapping database information to R. 
#' 
#' @param series character name key of 
#' @param map_list list to represent key value mapping. Could also be of class miro. 
#' @param meta_env an environment that already holds meta information and should be extended. 
#' Defaults to NULL in which case it creates and returns a new environment.
#' @param overwrite_objects logical should the entire existing meta information be overwritten inside
#' the environment? Defaults to FALSE
#' @param overwrite_elements logical should single matching elements of a meta information objectes be overwritten. 
#' Defaults to TRUE. 
#' @export
addMetaInformation <- function(series,map_list,
                   meta_env = NULL,
                   overwrite_objects = F,
                   overwrite_elements = T){
  # sanity check
  stopifnot(is.list(map_list))
  # check if environment exists, 
  # if not create it and put 
  # the meta information in there, stored
  # under the series' name.
  
  # general adjustment, add class
  # meta informaiton for R objects.
  class(map_list) <- c('miro','list')
  
  # remove empty elements from a list 
  # this can be important for generically
  # created meta information
  map_list[map_list == ''] <- NULL
  if(length(map_list) == 0) map_list <- NULL
    
  if(is.null(meta_env)){
    meta_env <- new.env()
    if(!is.null(map_list)){
      meta_env[[series]] <- map_list
    }    
  } else {
    # if environment exists we need to check
    # whether the object exists and if so
    # whether it needs be overwritten
    if(overwrite_objects){
      if(!is.null(map_list)){
        meta_env[[series]] <- map_list
      }
    } else {
      if(!is.null(meta_env[[series]])){
        
        # elements that are not in the current list should be 
        # attached
        elements_in_old <- (names(map_list) %in% names(meta_env[[series]]))
        new_elements <- map_list[!elements_in_old]
        meta_env[[series]] <- c(meta_env[[series]],new_elements)
        
        if(overwrite_elements & length(map_list[elements_in_old]) != 0){
          meta_env[[series]][names(map_list[elements_in_old])] <- map_list[elements_in_old]    
        }
      } else {
        meta_env[[series]] <- map_list
      }
    }
  }
  
  class(meta_env) <- c('meta_env','environment')
  meta_env
}

