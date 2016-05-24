#' Creates a subsetted observation table 
#'
#' Subsets the obs.table dataframe supplied to only contain the observations of
#' interest.
#'  
#' @param obs.table dataframe of observation records with fields object,
#'   Region.Label, and Sample.Label which give links to sample.table,
#'   region.table and the data records used in \code{model}
#' @param ddf.data dataframe containing the observations
#' @param subset.variable variable name supplied as a character
#' @param subset.value character value on which to subset the data
#' @return dataframe containing the subset of the obs.table
#' @note Internal function not intended to be called by user.
#' @author Laura Marshall
#' @keywords utility
#'
create.obs.table <- function(obs.table, ddf.data, subset.variable, subset.value){

  object.index <- which(ddf.data[[subset.variable]] == subset.value)
  object.id <- ddf.data[object.index,]$object
  
  obs.table.subset <- obs.table[which(obs.table$object%in%object.id),]  
  return(obs.table.subset)
}