#' Renumbers the object IDs for the duplicate observations generated when
#' bootstrapping
#'
#' Find the largest object ID and renumbers all duplicate IDs starting form 
#' this value. The information for the duplicates is also added to the obs.table
#'  
#' @param ddf.dat dataframe containing a single dataset with duplicate 
#'   observations
#' @param obs.table dataframe of observation records with fields object,
#'   Region.Label, and Sample.Label which give links to sample.table,
#'   region.table and the data records used in \code{model}
#' @param double.observer boolean indicating if it is a double observer survey
#' @return list with 2 elements:
#'   ddf.dat dataframe containing a single dataset with new and unique 
#'     observation IDs 
#'   obs.table the updated obs.table dataframe containing the new
#'     observation IDs
#' @note Internal function not intended to be called by user.
#' @author Laura Marshall
#' @seealso \code{resample.data}
#' @keywords data manipulation
#'
renumber.duplicates <- function(ddf.dat, obs.table, double.observer){
# renumber.duplicates function to renumbers the object IDs for the duplicate observations generated when
# bootstrapping.
#
# Arguments:
#   ddf.dat dataframe containing a single dataset with duplicate IDs 
#   obs.table dataframe of observation records
# 
# Value: list with 2 elements
#   ddf.dat updated dataframe containing a single dataset  
#   obs.table updated obs.table dataframe
#
# Function Call: none
#
  #find largest object ID
  next.id <- max(unique(obs.table$object)) + 1
  
  #get table of object IDs
  object.table <- table(ddf.dat$object)
  
  if(!double.observer){
    #get names of duplicated sightings
    duplicate.names <- names(which(object.table > 1))
    
    #for each duplicate name
    for(obj in seq(along=duplicate.names)){
      #find indexes of duplicates
      duplicate.index <- which(ddf.dat$object == duplicate.names[obj])
      #for each duplicate except the first
      for(dup in seq(along = duplicate.index)[-1]){
        #renumber object in ddf.dat
        ddf.dat$object[duplicate.index[dup]] <- next.id      
        #add row to the obs.table
        obs.table <- rbind(obs.table, obs.table[obs.table$object == duplicate.names[obj],])
        obs.table$object[nrow(obs.table)] <- next.id
        #increment counter
        next.id <- next.id + 1 
      } 
    } 
  }else{
    #get names of duplicated sightings
    duplicate.names <- names(which(object.table > 2))
    
    #for each duplicate name
    for(obj in seq(along=duplicate.names)){
      #find indexes of duplicates
      duplicate.index <- which(ddf.dat$object == duplicate.names[obj])
      #for each duplicate except the first 2
      dup <- 3
      while(dup < length(duplicate.index)){
        #renumber objects in ddf.dat
        ddf.dat$object[duplicate.index[dup]] <- next.id
        dup <- dup + 1 
        ddf.dat$object[duplicate.index[dup]] <- next.id   
        dup <- dup + 1  
        #add row to the obs.table
        obs.table <- rbind(obs.table, obs.table[obs.table$object == duplicate.names[obj],])
        obs.table$object[nrow(obs.table)] <- next.id
        #increment counter
        next.id <- next.id + 1 
      } 
    }
  }      
  return(list(ddf.dat = ddf.dat, obs.table = obs.table))
}