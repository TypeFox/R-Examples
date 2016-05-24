#' Function to return geochronological data from records.
#'
#' Using the dataset ID, return all geochronological data associated with the dataID.  At present,
#'    only returns the dataset in an unparsed format, not as a data table.   This function will only download one dataset at a time.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr GET content
#' @param x A numeric dataset ID or a vector of numeric dataset IDs, or an object of class of class \code{site}, \code{dataset}, \code{dataset_list}, \code{download} or \code{download_list} for which geochrons are required.
#' @param verbose logical; should messages on API call be printed?
#' 
#' @author Simon J. Goring \email{simon.j.goring@@gmail.com}
#' @return This command returns either an object of class \code{"try-error"}' (see \code{\link{try}}) definined by the error returned
#'    from the Neotoma API call, or a \code{geochronologic} object, which is a list with two components, a \code{dataset} and a geochronology table, a \code{data.frame} with the following components:
#'
#'  \item{ \code{sample.id} }{A unique identifier for the geochronological unit.}
#'  \item{ \code{age.type} }{String.  The age type, one of calendar years, radiocarbon years, etc.}
#'  \item{ \code{age} }{Dated age of the material.}
#'  \item{ \code{e.older} }{The older error limit of the age value.  Commonly 1 standard deviation.}
#'  \item{ \code{e.young} }{The younger error limit of the age value.}
#'  \item{ \code{delta13C} }{The measured or assumed delta13C value for radiocarbon dates, if provided.}
#'  \item{ \code{material.dated} }{A table describing the collection, including dataset information, PI data compatable with \code{\link{get_contact}} and site data compatable with \code{\link{get_site}}.}
#'  \item{ \code{geo.chron.type} }{Text string, type of geochronological analysis, i.e., Radiocarbon dating, luminesence.}
#'  \item{ \code{notes} }{Text string}
#'  \item{ \code{infinite} }{Boolean, does the dated material return an "infinte" date?}
#'
#'  A full data object containing all the relevant geochronological data available for a dataset.
#' @examples \dontrun{
#' #  Search for sites with "Pseudotsuga" pollen that are older than 8kyr BP and
#' #  find the relevant radiocarbon ages associated with the cores.
#' #  Are some time periods better dated than others?
#' t8kyr.datasets <- get_dataset(taxonname='*Pseudotsuga*', loc=c(-150, 20, -100, 60),
#'                               ageyoung = 8000)
#'
#' #  Returns 116 records (as of 13/07/2015).  These are the pollen records though, we want the sites:
#' geochron.records <- get_geochron(get_site(t8kyr.datasets))
#'
#' #  We want to extract all the radiocarbon ages from the records:
#'
#' get_ages <- function(x){
#'   any.ages <- try(x[[2]]$age[x[[2]]$age.type == 'Radiocarbon years BP'])
#'   if(class(any.ages) == 'try-error') output <- NA
#'   if(!class(any.ages) == 'try-error') output <- unlist(any.ages)
#'   output
#' }
#'
#' radio.chron <- unlist(sapply(geochron.records, get_ages))
#'
#' hist(radio.chron[radio.chron<40000], breaks=seq(0, 40000, by = 500),
#'      main = 'Distribution of radiocarbon dates for Pseudotsuga records',
#'      xlab = 'Radiocarbon date (14C years before 1950)')
#' }
#'
#' @references
#' Neotoma Project Website: http://www.neotomadb.org
#' API Reference:  http://api.neotomadb.org/doc/resources/contacts
#' @keywords IO connection
#' @export
get_geochron <- function(x, verbose = TRUE){
  UseMethod('get_geochron')
}

#' @importFrom jsonlite fromJSON
#' @importFrom httr GET content
#' @export
get_geochron.default <- function(x, verbose = TRUE){

  #  If it doesn't get passed through the other methods x needs to be numeric.
  if (!is.numeric(x)){
    stop('datasetid must be numeric.')
  }

  # Get sample is a function because we can now get
  # one or more geochronologies at a time.
  get_sample <- function(x){
    
    base.uri <- 'http://api.neotomadb.org/v1/apps/geochronologies'
    
    # query Neotoma for data set
    neotoma_content <- httr::content(httr::GET(paste0(base.uri, '/?datasetid=', x)), as = "text")
    if (identical(neotoma_content, "")) stop("")
    aa <- jsonlite::fromJSON(neotoma_content, simplifyVector = FALSE)
    
    # Might as well check here for error and bail
    if (inherits(aa, "try-error"))
        return(aa)

    # if no error continue processing
    if (isTRUE(all.equal(aa[[1]], 0))) {
      # The API did not return a record, or returned an error.
        stop(paste('Server returned an error message:\n', aa[[2]]),
             call. = FALSE)
    }

    if (isTRUE(all.equal(aa[[1]], 1) & length(aa[[2]]$Samples) == 0)) {
      # The API returned a record, but the record did not
      # have associated geochronology information.
      stop(paste0('No geochronological record is associated with dataset ',x,'.'),
           call. = FALSE)
    }

    if (isTRUE(all.equal(aa[[1]], 1) & length(aa[[2]]) > 0)) {
      # The API returned a record with geochron data.
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
      
      aa <- lapply(aa, function(x)rep_NULL(x))
      
      if (verbose) {
          message(strwrap(paste0("API call was successful.")))
      }

      # If there are actual stratigraphic samples
      # with data in the dataset returned.
      
      # We have to pull the dataset information from the `download`:
      dl <- try(jsonlite::fromJSON(paste0('http://api.neotomadb.org/v1/data/downloads/', x)))[[2]]
      
      dl <- lapply(dl, function(x)rep_NULL(x))
      
      dataset <- list(
        site.data = data.frame(site.id = dl$Site$SiteID,
                               site.name = dl$Site$SiteName,
                               long = mean(unlist(dl$Site[c('LongitudeWest', 'LongitudeEast')]),
                                           na.rm = TRUE),
                               lat = mean(unlist(dl$Site[c('LatitudeNorth', 'LatitudeSouth')]),
                                          na.rm = TRUE),
                               elev = dl$Site$Altitude,
                               description = dl$Site$SiteDescription,
                               long.acc = abs(dl$Site$LongitudeWest - dl$Site$LongitudeEast),
                               lat.acc = abs(dl$Site$LatitudeNorth - dl$Site$LatitudeSouth),
                               row.names = dl$Site$SiteName,
                               stringsAsFactors = FALSE),
        dataset.meta = data.frame(dataset.id = dl$DatasetID,
                                  dataset.name = dl$DatasetName,
                                  collection.type = dl$CollUnitType,
                                  collection.handle = dl$CollUnitHandle,
                                  dataset.type =  dl$DatasetType,
                                  stringsAsFactors = FALSE),
        pi.data = do.call(rbind.data.frame,
                          dl$DatasetPIs),
        submission = data.frame(submission.date = strptime(dl$NeotomaLastSub,
                                                           '%m/%d/%Y'),
                                submission.type = 'Last submission to Neotoma',
                                stringsAsFactors=FALSE),
        access.date = Sys.time())
      
      class(dataset) <- c('dataset', 'list')
      
      pull.rec <- function(x){

        data.frame(sample.id = x$SampleID,
                   depth   = x$Depth,
                   thickness = x$Thickness,
                   age.type = x$AgeType,
                   age = x$Age,
                   e.older = x$ErrorOlder,
                   e.young = x$ErrorYounger,
                   delta13C = x$Delta13C,
                   lab.no = x$LabNumber,
                   material.dated = x$MaterialDated,
                   geo.chron.type = x$GeochronType,
                   notes = x$Notes,
                   infinite = x$Infinite,
                   stringsAsFactors = FALSE)
      }

      out <- list(dataset = dataset, 
                  geochron   = do.call(rbind.data.frame, lapply(aa[[1]], pull.rec)))
      class(out) <- c('geochronologic', 'list')
      
      out
    }

    out
  }
  
  out <- lapply(x, function(x)try(get_sample(x)))
  for(i in length(out):1){if('try-error' %in% class(out[[i]]))out[[i]] <- NULL}
  
  if(length(out) == 0){
    # It's possible that we had some successes, if not then we need to return
    # an error.
    stop('There were no geochronological records returned.', call.=FALSE)
  }
  
  class(out) <- c('geochronologic_list', 'list')
  
  out

}

#' @export
get_geochron.dataset <- function(x, verbose = TRUE){
  
  # Updated the processing here. There is no need to be fiddling with
  # call. Use missing() to check for presence of argument
  # and then process as per usual
  
  datasetid <- x$dataset.meta$dataset.id
  
  if(!x$dataset.meta$dataset.type %in% 'geochronologic'){
    stop(paste0('The dataset ID ', x$dataset.meta$dataset.id,
                   ' is not associated with a geochronology object.'))
  } else {
    geochron <- get_geochron(datasetid)[[1]]
  }
  
  geochron[[1]] <- x
  
  class(geochron) <- c('geochronologic', 'list')
  
  geochron
  
}

#' @export
get_geochron.dataset_list <- function(x, verbose = TRUE){
  
  # Updated the processing here. There is no need to be fiddling with
  # call. Use missing() to check for presence of argument
  # and then process as per usual
  
  dataset.types <- unlist(lapply(x, FUN=function(x)x$dataset$dataset.type))
  
  if(any(!dataset.types%in%'geochronologic')){
    if(all(!dataset.types%in%'geochronologic')){
      stop('This set contains no geochronological datasets.  Use get_download instead.')
    } else {
      message('This dataset contains records that are not geochronological datasets.  Only geochronological datasets will be returned.')
      x <- x[dataset.types %in% 'geochronologic']
      
      if(length(x) > 1){
        class(x) <- c('dataset_list', 'list')
        
        aa <- lapply(x, function(y){
          out <- get_geochron(y)
          out
        })
        
        class(aa) <- c('geochronologic_list', 'list')
        
      } else {
        x <- x[[1]]
        aa <- get_geochron(x)
      }
    }
  }
  
  aa
}

#' @export
get_geochron.site <- function(x, verbose = TRUE){
  
  dataset <- get_dataset(x)
  
  dataset.types <- unlist(lapply(dataset, FUN=function(x)x$dataset.meta$dataset.type))
  
  if(any(!dataset.types%in%'geochronologic')){
    if(all(!dataset.types%in%'geochronologic')){
      stop('This set contains no geochronological datasets.  Use get_download instead.')
    } else {
      message('This dataset contains records that are not geochronological datasets.  Only geochronological datasets will be returned.')
      dataset <- dataset[dataset.types %in% 'geochronologic']
      class(dataset) <- c('dataset_list', 'list')
    }
  }
  
  aa <- lapply(dataset, function(y){
    out <- get_geochron(y)
    out
  })
  
  class(aa) <- c('geochronologic_list', 'list')
  
  aa

}

