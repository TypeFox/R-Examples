#' @title Obtain dataset information from the Neotoma Paleoecological Database or an existing object.
#' @description A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr content GET
#' @param x An optional value, either a \code{numeric} site ID or object of class \code{download}, \code{download_list} or \code{site}.
#' @param datasettype A character string corresponding to one of the allowed dataset types in the Neotoma Database.  Allowed types include: \code{"geochronologic"}, \code{"loss-on-ignition"}, \code{"pollen"}, \code{"plant macrofossils"}, \code{"vertebrate fauna"}, \code{"mollusks"}, and \code{"pollen surface sample"}.
#' @param piid Numeric value for the Principle Investigator's ID number.
#' @param altmin Numeric value indicating the minimum altitude for the site (can be used alone or with \code{altmax}).
#' @param altmax Numeric value indicating the maximum altitude for the site (can be used alone or with \code{altmin}).
#' @param loc A numeric vector \code{c(lonW, latS, lonE, latN)} representing the bounding box within which to search for sites.  The convention here is to use negative values for longitudes west of Grewnwich or longitudes south of the equator
#' @param gpid A character string or numeric value, must correspond to a valid geopolitical identity in the Neotoma Database.  Use get.tables('GeoPoliticalUnits') for a list of acceptable values, or link here: \url{http://api.neotomadb.org/apdx/geopol.htm}
#' @param taxonids A numeric identifier for the taxon.  See \code{\link{get_table}} and use \code{get_tables('Taxa')} for a list of acceptable values.
#' @param taxonname A character string corresponding to a valid taxon identity in the Neotoma Database.  See \code{\link{get_table}} and use \code{get_table('Taxa')} for a list of acceptable values.
#' @param ageold The oldest date acceptable for the search (in years before present).
#' @param ageyoung The youngest date acceptable for the search.
#' @param ageof If a taxon ID or taxon name is defined this parameter must be set to \code{"taxon"}, otherwise it may refer to \code{"sample"}, in which case the age bounds are for any samples within datasets or \code{"dataset"} if you want only datasets that are within the bounds of ageold and ageyoung.
#' @param subdate Date of dataset submission, either YYYY-MM-DD or MM-DD-YYYY.
#'
#' @author Simon J. Goring \email{simon.j.goring@@gmail.com}
#' @return More details on the use of these parameters can be obtained from
#'    \url{http://api.neotomadb.org/doc/resources/datasets}.
#'
#'    A list of class `dataset_list`, with each item corresponding to an individual record.
#'    Searches that return no items will result in a NULL value being returned.
#'    Otherwise each list item (each dataset record) includes the following components:
#'
#'  \item{ \code{dataset.id} }{Unique database record identifier for the dataset.}
#'  \item{ \code{dataset.name}  }{Name of the dataset; not commonly used.}
#'  \item{ \code{CollUnitHandle}  }{Code name of the Collection Unit with which the dataset is associated. This code may be up to 10 characters. Data are frequently distributed by Collection Unit, and the Handle is used for file names.}
#'  \item{ \code{CollUnitID}  }{Unique database record identifier for the collection unit.}
#'  \item{ \code{CollType}  }{The collection type. Types include cores, sections, excavations, and animal middens.}
#'  \item{ \code{DatasetType}  }{The dataset type, such as: geochronologic, loss-on-ignition, pollen, plant macrofossils, vertebrate fauna, etc.}
#'  \item{ \code{AgeOldest}  }{The oldest of all sample ages (in calendar years before present) in the dataset.}
#'  \item{ \code{AgeYoungest}  }{The youngest of all sample ages (in calendar years before present) in the dataset.}
#'  \item{ \code{SubDates}  }{An array of objects that describe dataset submission events.  If multiple submissions occured then this is a table.}
#'  \item{ \code{DatasetPIs}  }{An array of objects that describe Principal Investigators associated with a dataset.}
#'  \item{ \code{Site}  }{An object describing the site where the dataset samples were taken.}
#'
#' @examples \dontrun{
#' # Search for sites with "Thuja" pollen that are older than 8kyr BP and
#' # that are on the west coast of North America:
#' t8kyr.datasets <- get_dataset(taxonname='Thuja*', loc=c(-150, 20, -100, 60), ageyoung = 8000)
#'
#' # Search for vertebrate fossils in Canada (gpid: 756) within the last 2kyr.
#' gpids <- get_table(table.name='GeoPoliticalUnits')
#' canID <- gpids[which(gpids$GeoPoliticalName == 'Canada'),1]
#'
#' v2kyr.datasets <- get_dataset(datasettype='vertebrate fauna', gpid=canID, ageold = 2000)
#' }
#' @references
#' Neotoma Project Website: http://www.neotomadb.org
#' API Reference:  http://api.neotomadb.org/doc/resources/contacts
#' @keywords IO connection
#' @export
#'
get_dataset <- function(x, datasettype, piid, altmin, altmax, loc, gpid, taxonids, taxonname, ageold, ageyoung, ageof, subdate) {
  UseMethod('get_dataset')
}

#' @title Obtain dataset information from the Neotoma Paleoecological Database or an existing object.
#' @description A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr GET content
#' @param x A numeric value corresponding to the site ID.
#' @param datasettype A character string corresponding to one of the allowed dataset types in the Neotoma Database.  Allowed types include: \code{"geochronologic"}, \code{"loss-on-ignition"}, \code{"pollen"}, \code{"plant macrofossils"}, \code{"vertebrate fauna"}, \code{"mollusks"}, and \code{"pollen surface sample"}.
#' @param piid Numeric value for the Principle Investigator's ID number.
#' @param altmin Numeric value indicating the minimum altitude for the site (can be used alone or with \code{altmax}).
#' @param altmax Numeric value indicating the maximum altitude for the site (can be used alone or with \code{altmin}).
#' @param loc A numeric vector \code{c(lonW, latS, lonE, latN)} representing the bounding box within which to search for sites.  The convention here is to use negative values for longitudes west of Grewnwich or longitudes south of the equator
#' @param gpid A character string or numeric value, must correspond to a valid geopolitical identity in the Neotoma Database.  Use get.tables('GeoPoliticalUnits') for a list of acceptable values, or link here: \url{http://api.neotomadb.org/apdx/geopol.htm}
#' @param taxonids A numeric identifier for the taxon.  See \code{\link{get_table}} and use \code{get_tables('Taxa')} for a list of acceptable values.
#' @param taxonname A character string corresponding to a valid taxon identity in the Neotoma Database.  See \code{\link{get_table}} and use \code{get_table('Taxa')} for a list of acceptable values.
#' @param ageold The oldest date acceptable for the search (in years before present).
#' @param ageyoung The youngest date acceptable for the search.
#' @param ageof If a taxon ID or taxon name is defined this parameter must be set to \code{"taxon"}, otherwise it may refer to \code{"sample"}, in which case the age bounds are for any samples within datasets or \code{"dataset"} if you want only datasets that are within the bounds of ageold and ageyoung.
#' @param subdate Date of dataset submission, either YYYY-MM-DD or MM-DD-YYYY.
#' @export
get_dataset.default <- function(x, datasettype, piid, altmin, altmax, loc, gpid, taxonids, taxonname, ageold, ageyoung, ageof, subdate) {
  # The issue here is that these objects
  # have multiple tables of multiple lengths.

  base.uri <- 'http://api.neotomadb.org/v1/data/datasets'

  cl <- as.list(match.call())
  cl[[1]] <- NULL
  cl <- lapply(cl, eval, envir = parent.frame())
  
  if ('x' %in% names(cl)) {
    names(cl)[which(names(cl) == 'x')] <- 'siteid'
  }

  #  Pass the parameters to param_check to make sure everything is kosher.
  error_test <- param_check(cl)
  if (error_test[[2]]$flag == 1) {
    stop(paste0(unlist(error_test[[2]]$message), collapse='\n  '))
  } else {
    cl <- error_test[[1]]
  }

  neotoma_content <- httr::content(httr::GET(base.uri, query = cl), as = "text")
  if (identical(neotoma_content, "")) stop("")
  aa <- jsonlite::fromJSON(neotoma_content, simplifyVector = FALSE)
  
  if (aa[[1]] == 0) {
    stop(paste('Server returned an error message:\n', aa[[2]]), call. = FALSE)
  }
  if (aa[[1]] == 1) {
    output <- aa[[2]]
    
    rep_NULL <- function(x) { 
      if (is.null(x)) {NA}
      else{
        if (class(x) == 'list') {
          lapply(x, rep_NULL)
        } else {
          return(x)
        }
      }
    }
    
    output <- lapply(output, function(x)rep_NULL(x))
    
    if (length(output) == 0) {
      warning('The criteria used returned 0 sample sites. Returning NULL.')
      return(NULL)
    }
    
    if (length(aa[[2]]) > 1) {
      message(paste('The API call was successful, you have returned ',
                    length(output), ' records.\n', sep = ''))
    } else {
      message(paste('The API call was successful, you have returned ',
                    length(output), ' record.\n', sep = ''))
    }
  }


  if (inherits(output, "try-error")) {
      new.output <- output
  } else {
      new.output <- lapply(output, function(x) {
          new.output <- list()
          new.output$site.data <- data.frame(site.id = x$Site$SiteID,
                                             site.name = x$Site$SiteName,
                                             long = mean(unlist(x$Site[c('LongitudeWest', 'LongitudeEast')]),
                                             na.rm = TRUE),
                                             lat = mean(unlist(x$Site[c('LatitudeNorth', 'LatitudeSouth')]),
                                             na.rm = TRUE),
                                             elev = x$Site$Altitude,
                                             description = x$Site$SiteDescription,
                                             long.acc = abs(x$Site$LongitudeWest - x$Site$LongitudeEast),
                                             lat.acc = abs(x$Site$LatitudeNorth - x$Site$LatitudeSouth),
                                             row.names = x$Site$SiteName,
                                             stringsAsFactors = FALSE)

          class(new.output$site.data) <- c('site', 'data.frame')

          if ('CollType' %in% names(x)) {x$CollUnitType <- x$CollType} # This is a fix for a very specific issue we were having.
          
          new.output$dataset.meta <- data.frame(dataset.id = ifelse(class(x$DatasetID) == 'logical',
                                                                    NA, x$DatasetID),
                                                dataset.name = ifelse(class(x$DatasetName) == 'logical',
                                                                      NA, x$DatasetName),
                                                collection.type = ifelse(class(x$CollUnitType) == 'logical',
                                                                         NA, x$CollUnitType),
                                                collection.handle = ifelse(class(x$CollUnitHandle) == 'logical',
                                                                           NA, x$CollUnitHandle),
                                                dataset.type = ifelse(class(x$DatasetType) == 'logical',
                                                                      NA, x$DatasetType),
                                                stringsAsFactors = FALSE)
          if (class(x$DatasetPIs) == 'logical') { 
            new.output$pi.data <- NA
          } else {
            new.output$pi.data <- do.call(rbind.data.frame, x$DatasetPIs)
            rownames(new.output$pi.data) <- NULL
          }

          sub.test <- try(do.call(rbind.data.frame, x$SubDates))

          if (length(sub.test) > 0) {
              colnames(sub.test) <- c("submission.date",  "submission.type")
              sub.test$submission.date <- as.character(sub.test$submission.date)
              sub.test$submission.type <- as.character(sub.test$submission.type)
          }

          new.output$submission <- sub.test

          new.output$access.date = Sys.time()

          class(new.output) <- c('dataset', 'list')
          new.output})

  }
  
  names(new.output) <- sapply(lapply(new.output, '[[', 'dataset.meta'), 
                              '[[', 'dataset.id')
  
  class(new.output) <- c('dataset_list', 'list')

  new.output

}

#' @title Obtain dataset information from an existing \code{site} object.
#' @description A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
#'
#' @param x An object of class \code{site}.
#' @param ... objects passed from the generic.  Not used in the call.
#' @importFrom jsonlite fromJSON
#' @importFrom httr GET content
#' @export
get_dataset.site <- function(x, ...) {

  rep_NULL <- function(x) { 
    if (is.null(x) | length(x) == 0) {NA}
    else{
      if (class(x) == 'list') {
        lapply(x, rep_NULL)
      } else {
        return(x)
      }
    }
  }
  
  
  pull_site <- function(siteid) {
    
    base.uri <- 'http://api.neotomadb.org/v1/data/datasets/?siteid='
    
    neotoma_content <- httr::content(httr::GET(paste0(base.uri, siteid)), as = "text")
    if (identical(neotoma_content, "")) stop("")
    aa <- jsonlite::fromJSON(neotoma_content, simplifyVector = FALSE)

    if (aa[[1]] == 0) {
      stop(paste('Server returned an error message:\n', aa[[2]]), call. = FALSE)
    }
    if (aa[[1]] == 1) {
      output <- aa[[2]]
      # replace NULL values:
      
    }

    if (length(output) == 0) {
      warning('The criteria used returned 0 sample sites. Returning NULL.')
      return(NULL)
    }
    
    new.output <- lapply(output, function(x) {
      new.output <- list()
      
      x <- rep_NULL(x)
      
      new.output$site.data <- data.frame(site.id = x$Site$SiteID,
                                         site.name = x$Site$SiteName,
                                         long = mean(unlist(x$Site[c('LongitudeWest', 'LongitudeEast')]),
                                                     na.rm = TRUE),
                                         lat = mean(unlist(x$Site[c('LatitudeNorth', 'LatitudeSouth')]),
                                                    na.rm = TRUE),
                                         elev = x$Site$Altitude,
                                         description = x$Site$SiteDescription,
                                         long.acc = abs(x$Site$LongitudeWest - x$Site$LongitudeEast),
                                         lat.acc = abs(x$Site$LatitudeNorth - x$Site$LatitudeSouth),
                                         row.names = x$Site$SiteName,
                                         stringsAsFactors = FALSE)
      
      new.output$dataset.meta <- data.frame(dataset.id = ifelse(class(x$DatasetID) == 'logical',
                                                                NA, x$DatasetID),
                                            dataset.name = ifelse(class(x$DatasetName) == 'logical',
                                                                  NA, x$DatasetName),
                                            collection.type = ifelse(class(x$CollUnitType) == 'logical',
                                                                     NA, x$CollUnitType),
                                            collection.handle = ifelse(class(x$CollUnitHandle) == 'logical',
                                                                       NA, x$CollUnitHandle),
                                            dataset.type = ifelse(class(x$DatasetType) == 'logical',
                                                                  NA, x$DatasetType),
                                            stringsAsFactors = FALSE)
      if (class(x$DatasetPIs) == 'logical') { 
        new.output$pi.data <- NA
      } else {
        new.output$pi.data <- do.call(rbind.data.frame, x$DatasetPIs)
        rownames(new.output$pi.data) <- NULL
      }
      
      sub.test <- try(do.call(rbind.data.frame, x$SubDates), silent=TRUE)

      if (length(sub.test) > 0 & !"try-error" %in% class(sub.test)) {
        colnames(sub.test) <- c("SubmissionDate",  "SubmissionType")
      } else {
        sub.test <- data.frame(SubmissionDate = NA, SubmissionType = NA)
      }

      new.output$submission <- sub.test

      new.output$access.date = Sys.time()

      class(new.output) <- c('dataset', 'list')
      
      new.output})

    new.output
  }

  new.output <- unlist(lapply(x$site.id,pull_site), recursive=FALSE)
  
  class(new.output) <- c('dataset_list', 'list')
  
  new.output

}

#' @title Obtain dataset information from an existing \code{download} object.
#' @description A function to access a \code{dataset} within a \code{download} object.
#'
#' @param x An object of class \code{download}.
#' @param ... objects passed from the generic.  Not used in the call.
#' @export
get_dataset.download <- function(x, ...) {
  # Just pull the dataset out of the download.
  output <- list(x$dataset)

  names(output) <- output[[1]]$dataset.meta$dataset.id

  class(output[[1]]) <- c('dataset', 'list')

  class(output) <- c('dataset_list', 'list')
  return(output)
}

#' @title Obtain dataset information from a \code{download_list}.
#' @description A function to return datasets corresponding to the objects within a \code{download_list}.
#'
#' @param x An object of class \code{download_list}.
#' @param ... objects passed from the generic.  Not used in the call.
#' @export
get_dataset.download_list <- function(x, ...) {

  # Just pull the dataset out of the download and reassign classes:
  output <- lapply(x, FUN=function(y) {
    dataset <- y$dataset
    class(dataset) <- c('dataset', 'list')
    dataset })

  names(output) <- sapply(lapply(output, '[[', 'dataset.meta'), '[[', 'dataset.id')
  
  class(output) <- c('dataset_list', 'list')

  output
}

#' @title Obtain dataset information from an object of class \code{geochronologic}.
#' @description A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
#'
#' @param x An object of class \code{geochronologic}.
#' @param ... objects passed from the generic.  Not used in the call.
#' @export
get_dataset.geochronologic <- function(x, ...) {
  x[[1]]
}

#' @title Obtain dataset information from an object of class \code{geochronologic_list}.
#' @description A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
#'
#' @param x An object of class \code{geochronologic_list}.
#' @param ... objects passed from the generic.  Not used in the call.
#' @export
get_dataset.geochronologic_list <- function(x, ...) {
  out <- lapply(x, function(y)y[[1]])
  class(out) <- c('dataset_list', 'list')
  out
}
