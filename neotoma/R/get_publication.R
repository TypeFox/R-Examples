#' @title A function to get publications for sites or datasets in the Neotoma Database using the API.
#'
#' @description The function takes the parameters, defined by the user, and returns a table with publication information from the Neotoma Paleoecological Database.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr GET content
#' @param x Numeric Publication ID value, either from \code{\link{get_dataset}} or known.
#' @param contactid Numeric Contact ID value, either from \code{\link{get_dataset}} or \code{\link{get_contact}}
#' @param datasetid Numeric Dataset ID, known or from \code{\link{get_dataset}}
#' @param author Character string for full or partial author's name.  Can include wildcards such as 'Smit*' for all names beginning with 'Smit'.
#' @param pubtype Character string, one of eleven allowable types, see \code{\link{get_table}}. For a list of allowed types run \code{get_table("PublicationTypes")}.
#' @param year Numeric publication year.
#' @param search A character string to search for within the article citation.
#'
#' @author Simon J. Goring \email{simon.j.goring@@gmail.com}
#' @return A list is returned with two data frame components:
#'
#'  \item{ \code{meta} }{A single row with Publication ID, type, year of publication and full citation.}
#'  \item{ \code{Authors} }{\code{data.frame} of author names, order and IDs, can be of variable length.}
#'
#' @examples \dontrun{
#' #  To find all publications from 1998:
#' year.cont <- get_publication(year = 1998)
#'
#' # To find all data contributors who have the last name "Smith"
#' smith.cont <- get_publication(author = 'Smith')
#' }
#' @references
#' Neotoma Project Website: http://www.neotomadb.org
#' API Reference:  http://api.neotomadb.org/doc/resources/contacts
#' @keywords IO connection
#' @export
#' 
get_publication<- function(x, contactid, datasetid, author,
                           pubtype, year, search){
  
  UseMethod('get_publication')

}


#' A function to get publications for sites or datasets in the Neotoma Database using the API.
#'
#' The function takes the parameters, defined by the user, and returns a table with
#'    publication information from the Neotoma Paleoecological Database.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr content GET
#' @param x Numeric Publication ID value, either from \code{\link{get_dataset}} or known.
#' @param contactid Numeric Contact ID value, either from \code{\link{get_dataset}} or \code{\link{get_contact}}
#' @param datasetid Numeric Dataset ID, known or from \code{\link{get_dataset}}
#' @param author Character string for full or partial author's name.  Can include wildcards such as 'Smit*' for all names beginning with 'Smit'.
#' @param pubtype Character string, one of eleven allowable types, see \code{\link{get_table}}. For a list of allowed types run \code{get_table("PublicationTypes")}.
#' @param year Numeric publication year.
#' @param search A character string to search for within the article citation.
#' @export
#' 
get_publication.default <- function(x, contactid, datasetid, author,
                            pubtype, year, search){

  base.uri <- 'http://api.neotomadb.org/v1/data/publications'

  cl <- as.list(match.call())
  cl[[1]] <- NULL
  
  if('x' %in% names(cl)){
    names(cl)[which(names(cl) == 'x')] <- 'pubid'
  }
  cl <- lapply(cl, eval, envir = parent.frame())

  #  Pass the parameters to param_check to make sure everything is kosher.
  error_test <- param_check(cl)
  if(error_test[[2]]$flag == 1){
    stop(paste0(unlist(error_test[[2]]$message), collapse='\n  '))
  } else {
    cl <- error_test[[1]]
  }
  
  
  neotoma_content <- httr::content(httr::GET(base.uri, query = cl), as = "text")
  if (identical(neotoma_content, "")) stop("")
  aa <- jsonlite::fromJSON(neotoma_content, simplifyVector = FALSE)
  
  if (aa[[1]] == 0){
    stop(paste('Server returned an error message:\n', aa[[2]]), call. = FALSE)
  }
  if (aa[[1]] == 1){
    aa <- aa[[2]]
    
    rep_NULL <- function(x){ 
      if(is.null(x)){NA}
      else{
        if(class(x) == 'list'){
          lapply(x, rep_NULL)
        } else {
          return(x)
        }
      }
    }
    
    # Clear NULLs from the output object & replace with NA values.
    aa <- lapply(aa, function(x)rep_NULL(x))
    
    if(length(aa) > 1){
      cat('The API call was successful, you have returned ',
          length(aa), 'records.\n')
    } else {
      cat('The API call was successful, you have returned  1 record.\n')
    }
  }

  if (class(aa) == 'try-error'){
    output <- NA
  } else {

    get_results <- function(x){
      
      if(length(x) == 0){
        output <- list(meta = data.frame(id = NA,
                                         pub.type = NA,
                                         year = NA,
                                         citation = NA,
                                         stringsAsFactors=FALSE))
        
        output$authors <- data.frame(contact.id = NA,
                                     order = NA,
                                     contact.name = NA,
                                               stringsAsFactors=FALSE)
        
      } else {
        
        output <- list(meta = data.frame(id = as.numeric(x$PublicationID),
                                         pub.type = x$PubType,
                                         year = as.numeric(x$Year),
                                         citation = x$Citation,
                                         stringsAsFactors=FALSE))
        
        output$authors <- do.call(rbind.data.frame,
          lapply(x$Authors, FUN=function(y){
          data.frame(ContactID = y$ContactID,
                     Order = y$Order,
                     ContactName = as.character(y$ContactName),
                     stringsAsFactors=FALSE)}))
      }
      
      output
    }

    if(length(aa) < 1) {
      output <- get_results(aa)
    } else{
      output <- lapply(aa, get_results)
    }
  }
  
  if('meta' %in% names(output)) output <- list(output)
  
  output
}

#' @title A function to get publications for datasets in the Neotoma Database using the API.
#' @description The function takes a \code{dataset} and returns a table with publication information from the Neotoma Paleoecological Database.
#'
#' @param x an object of class \code{dataset}.
#' @param ... objects passed from the generic.  Not used in the call.
#' @export
get_publication.dataset <- function(x, ... ){
  pubs <- get_publication(datasetid = x$dataset.meta$dataset.id)
  if('meta' %in% names(pubs)) pubs <- list(pubs)
  pubs
}

#' @title A function to get publications for dataset_lists in the Neotoma Database using the API.
#' @description The function takes a \code{dataset_list} and returns a table with publication information from the Neotoma Paleoecological Database.
#'
#' @param x an object of class \code{dataset_list}.
#' @param ... objects passed from the generic.  Not used in the call.
#' @export
get_publication.dataset_list <- function(x, ... ){
  ids <- sapply(x, function(y)y$dataset.meta$dataset.id)
  
  lapply(ids, function(x)get_publication(datasetid = x))
}

#' @title A function to get publications for downloads in the Neotoma Database using the API.
#' @description The function takes a \code{download} and returns a table with publication information from the Neotoma Paleoecological Database.
#'
#' @param x an object of class \code{download}.
#' @param ... objects passed from the generic.  Not used in the call.
#' @export
get_publication.download <- function(x, ... ){
  
  pubs <- get_publication(datasetid = x$dataset$dataset.meta$dataset.id)
  if('meta' %in% names(pubs)) pubs <- list(pubs)
  pubs
  
}

#' @title A function to get publications for datasets in the Neotoma Database using the API.
#' @description The function takes a \code{download_list} and returns a table with publication information from the Neotoma Paleoecological Database.
#'
#' @param x an object of class \code{download_list}.
#' @param ... objects passed from the generic.  Not used in the call.
#' @export
get_publication.download_list <- function(x, ... ){
  
  ids <- sapply(x, function(y)y$dataset$dataset.meta$dataset.id)
  
  lapply(ids, function(x)get_publication(datasetid = x))
  
}