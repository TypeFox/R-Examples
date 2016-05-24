#' @title Open a browser window to display a Neotoma dataset within the Neotoma Explorer
#' @description Using a \code{download} or \code{dataset} object, open up a browser window in the users default browser. Passing a \code{download_list} or \code{dataset_list} will open Neotoma Explorer with the first object and return a warning.
#'
#' @importFrom utils browseURL
#' @param x A \code{numeric} value for the dataset ID, a \code{dataset} or \code{download} object.
#'
#' @author Simon J. Goring \email{simon.j.goring@@gmail.com}
#' @return Returns a NULL value, opens a browser.
#' 
#' @examples \dontrun{
#' # Where are the XRF data?
#'
#' xrf.data <- get_dataset(datasettype='X-ray fluorescence (XRF)')
#' browse(xrf.data)
#'
#' }
#' @references
#' Neotoma Project Website: http://www.neotomadb.org
#' API Reference:  http://api.neotomadb.org/doc/resources/sites
#' @keywords IO connection
#' @export
browse <- function(x){
  UseMethod('browse', object = x)
}

#' @title Open a browser window to display a Neotoma dataset within the Neotoma Explorer
#' @description Using a \code{download} or \code{dataset} object, open up a browser window in the users default browser. Passing a \code{download_list} or \code{dataset_list} will open Neotoma Explorer with the first object and return a warning.
#'
#' @param x A numeric value with the dataset ID.
#' 
#' @export
browse.default <- function(x){
  if(length(x)>1){
    warning('Can only open one site at a time currently.  Opening the first dataset.')
  }
  utils::browseURL(paste0('http://apps.neotomadb.org/Explorer/?datasetid=', x))
  NULL
}


#' @title Open a browser window to display a Neotoma dataset within the Neotoma Explorer
#' @description Using a \code{download} or \code{dataset} object, open up a browser window in the users default browser. Passing a \code{download_list} or \code{dataset_list} will open Neotoma Explorer with the first object and return a warning.
#'
#' @param x A \code{dataset} object.
#' 
#' @export
browse.dataset <- function(x){
  browse(x$dataset.meta$dataset.id)
  NULL
}

#' @title Open a browser window to display a Neotoma dataset within the Neotoma Explorer
#' @description Using a \code{download} or \code{dataset} object, open up a browser window in the users default browser. Passing a \code{download_list} or \code{dataset_list} will open Neotoma Explorer with the first object and return a warning.
#'
#' @param x A \code{dataset_list} object.
#' 
#' @export
browse.dataset_list <- function(x){
  if(length(x) > 1){
    warning(paste0('This dataset_list has more than one site.  Currently the API only supports \n',
                   'displaying a single record at a time.  Displaying the first dataset in the list.'))
  }
  browse(x[[1]]$dataset.meta$dataset.id)
  NULL
}

#' @title Open a browser window to display a Neotoma dataset within the Neotoma Explorer
#' @description Using a \code{download} or \code{dataset} object, open up a browser window in the users default browser. Passing a \code{download_list} or \code{dataset_list} will open Neotoma Explorer with the first object and return a warning.
#'
#' @param x A \code{download} object.
#' 
#' @export
browse.download <- function(x){

  browse(x$dataset$dataset.meta$dataset.id)
  NULL
}

#' @title Open a browser window to display a Neotoma dataset within the Neotoma Explorer
#' @description Using a \code{download} or \code{dataset} object, open up a browser window in the users default browser. Passing a \code{download_list} or \code{dataset_list} will open Neotoma Explorer with the first object and return a warning.
#'
#' @param x A \code{download_list} object.
#' 
#' @export
browse.download_list <- function(x){
  if(length(x) > 1){
    warning(paste0('This download_list has more than one site.  Currently the API only supports \n',
                   'displaying a single record at a time.  Displaying the first download in the list.'))
  }
  browse(x[[1]]$dataset$dataset.meta$dataset.id)
  NULL
}