
#' Function to bind objects together into a longer object.
#'
#' From multiple \code{download*}s, \code{dataset*}s or \code{site}s, join them together into a single object.
#'
#' To support further synthesis and analysis \code{compile_download} works to transform a list
#' returned by \code{\link{get_download}} into a large data frame with columns for site and sample attributes
#' and also with the associated assemblage data at each sample depth.  This function also does the same for
#' single sites.
#'
#' @param x An object returned by one of the \code{get_*} commands for download, site or dataset.
#' @param  ... other objects of the same class.
#' @author Simon J. Goring \email{simon.j.goring@@gmail.com}
#' @return This command returns a larger list.
#'
#' @examples \dontrun{
#' #  Search for sites with "Thuja" pollen that are older than 8kyr BP and
#' #  that are on the west coast of North America:
#' t8kyr.poa <- get_dataset(taxonname='Thuja*', 
#'                          loc=c(-150, 20, -100, 60), ageyoung = 8000)
#' t8kyr.canis <- get_dataset(taxonname='Canis*', 
#'                            loc=c(-150, 20, -100, 60), ageyoung = 8000)
#'
#' t8kyr.co_site <- bind(t8kyr.poa, t8kyr.canis)
#' plot(t8kyr.co_site)
#' 
#' ####
#' # Download dataset data across four dataset types along the forest prairie boundary:
#' # We want to look at four different dataset types across a forest prairie boundary:
#' dataset_types <- c('ostracode surface sample',
#'                    'water chemistry',
#'                    'diatom surface sample',
#'                    'pollen surface sample')
#' 
#' # Run the 'get_dataset` function for each of the different dataset types 
#' dataset_lists <- lapply(dataset_types, 
#'                           function(x) { 
#'                             get_dataset(datasettype=x, 
#'                                         loc = c(-100,43,-92,48))
#'                                         })
#' 
#' # Using do.call here to make sure that I don't have to split the list out.
#' new_datasets <- do.call(bind, dataset_lists)
#' 
#' # And voila!
#' plot(new_datasets)
#' 
#' }
#' @references
#' Neotoma Project Website: http://www.neotomadb.org
#'
#' @keywords utilities
#' @export

bind <- function(x, ...) {
  
  inputs <- list(x, ...)
  
  if (!length(inputs)>1) {
    # If there's only one element, make sure it's a list and that the 
    # classes within the list are equivalent.
    inputs <- x
    if ('list' %in% class(inputs)) {
      list_class <- sapply(inputs, class)
      if (do.call(all.equal, as.list(list_class[1,]))) {
        # All datasets are of the same type:
        return(do.call(bind, inputs))
      } else {
        stop('All list elements must be of the same type.')
      }
    } else{
      stop('When passing a single object, that object must be a list of all the same data type.')
    }
  }
  
  classes <- sapply(inputs, function(x)unlist(class(x)))
  
  if (all(classes[1,] == 'download_list')) {
    new.list <- do.call(c, inputs)
    class(new.list) <- c('download_list', 'list')
    return(new.list)
  } else if (all(classes[1,] == 'dataset_list')) {
    new.list <- do.call(c, inputs)
    class(new.list) <- c('dataset_list', 'list')
    return(new.list)
  } else if (all(classes[1,] == 'site')) {
    new.list <- do.call(rbind.data.frame, inputs)
    class(new.list) <- c('site', 'data.frame')
    return(new.list)
  } else if (!all(classes[1,] == 'download_list') & all(classes[1,] %in% c('download_list', 'download'))) {
    #  Turn them into download_lists first.
    inputs <- lapply(inputs, function(x) {
      if (class(x)[1] == 'download') {
        x <- list(x)
        class(x) <- c('download_list', 'list')
      }
      x})
    
    new.list <- do.call(c, inputs)
    class(new.list) <- c('download_list', 'list')
    return(new.list)
    
  } else if (!all(classes[1,] == 'dataset_list') & all(classes[1,] %in% c('dataset_list', 'dataset'))) {
    #  Turn them into dataset_lists first.
    x <- lapply(inputs, function(x) {
      if (class(x)[1] == 'dataset') {
        x <- list(x)
        class(x) <- c('dataset_list', 'list')
      }
      x})

    new.list <- do.call(c, inputs)
    class(new.list) <- c('dataset_list', 'list')
    return(new.list)
    
  }
  else {
    stop('Objects must be of the same class (or related "_list" classes)')
  }
  
}