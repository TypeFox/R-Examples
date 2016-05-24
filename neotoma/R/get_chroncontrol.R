#' Function to return chronological control tables used to build age models.
#'
#' Using the dataset ID, return all records associated with the data.  At present,
#'    only returns the dataset in an unparsed format, not as a data table.   This function will only download one dataset at a time.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr GET content
#' @param x A single numeric chronology ID or a vector of numeric dataset IDs as returned by \code{\link{get_dataset}}.
#' @param verbose logical, should messages on API call be printed?
#' @param add logical, should this chron control be added to the download object?
#' @author Simon J. Goring \email{simon.j.goring@@gmail.com}
#' @return This command returns either an object of class  \code{"try-error"} containing the error returned
#'    from the Neotoma API call, or a full data object containing all the relevant information required to build either the default or prior chronology for a core.
#'    When \code{download} or \code{download_list} objects are passes, the user can \code{add} the chroncontrol to the
#'    \code{download} object explicitly, in which case the function will return a download with \code{chroncontrol} embedded.
#'    
#'    This is a list comprising the following items:
#'
#'  \item{ \code{chron.control} }{A table describing the collection, including dataset information, PI data compatable with \code{\link{get_contact}} and site data compatable with \code{\link{get_site}}.}
#'  \item{ \code{meta} }{Dataset information for the core, primarily the age-depth model and chronology.  In cases where multiple age models exist for a single record the most recent chronology is provided here.}
#'  
#'  If Neotoma returns empty content, either the control table or the associated metadata (which happens in approximately 25% of cases) then the data.frames are returned with NA content.
#'  
#' @examples \dontrun{
#' #  The point of pulling chronology tables is to re-build or examine the chronological
#' #  information that was used to build the age-depth model for the core.
#' 
#' }
#' @references
#' Neotoma Project Website: http://www.neotomadb.org
#' API Reference:  http://api.neotomadb.org/doc/resources/contacts
#' @keywords IO connection
#' @export
get_chroncontrol <- function(x, verbose = TRUE, add = FALSE) {
  UseMethod('get_chroncontrol')
}

#' @title Function to return chronological control tables from a chronologic ID.
#' @description Using the chronology ID, return the chron control table as a \code{data.frame}.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr GET content
#' @param x A single numeric chronology ID or a vector of numeric chronology IDs as returned by \code{get_datasets}.
#' @param add logical, should this chron control be added to the download object?
#' @param verbose logical; should messages on API call be printed?
#' @export
get_chroncontrol.default <- function(x, verbose = TRUE, add = FALSE) {
  
  # Updated the processing here. There is no need to be fiddling with
  # call. Use missing() to check for presence of argument
  # and then process as per usual
  base.uri <- 'http://api.neotomadb.org/v1/data/chronologies'
  
  if (missing(x)) {
    stop(paste(sQuote("chronologyid"), "(x) must be provided."))
  } else {
    if (!is.numeric(x))
      stop('chronology ID (x) must be numeric.')
  }
  
  # query Neotoma for data set
  neotoma_content <- httr::content(httr::GET(paste0(base.uri, '/', x)), as = "text")
  if (identical(neotoma_content, "")) stop("")
  aa <- jsonlite::fromJSON(neotoma_content, simplifyVector = FALSE)
  
  # Might as well check here for error and bail
  if (inherits(aa, "try-error"))
    return(aa)
  
  # if no error continue processing
  if (isTRUE(all.equal(aa[[1]], 0))) {
    stop(paste('Server returned an error message:\n', aa[[2]]),
         call. = FALSE)
  }
  
  if (isTRUE(all.equal(aa[[1]], 1))) {
    aa <- aa[[2]]
    
    rep_NULL <- function(x) { 
      if (is.null(x)) { NA }
      else{
        if (class(x) == 'list') {
          lapply(x, rep_NULL)
        } else {
          return(x)
        }
      }
    }
    
    aa <- lapply(aa, function(x)rep_NULL(x))
    
    if (verbose) {
      writeLines(strwrap(paste0("API call was successful.",
                                " Returned chronology.")))
    }
    
    # Here the goal is to reduce this list of lists to as
    # simple a set of matrices as possible.
    #  Some of the records do not contain a 'controls' table, so they only contain
    #  'meta':
    
    if (length(aa) == 0) {
      # This removes the need to make a more deeply nested if statement. . . 
      aa <- list(empty = 0) 
    }
    
    if ('controls' %in% names(aa[[1]])) {
      
      flattened <- lapply(lapply(aa, '[[', 'controls')[[1]], data.frame)
      control.table <- do.call(rbind.data.frame, flattened)
      
      control.table <- control.table[, c('Depth', 'Thickness',
                                         'Age', 'AgeYoungest', 'AgeOldest',
                                         'ControlType', 'ChronControlID')]
      
      colnames(control.table) <- c('depth', 'thickness', 'age', 
                                   'age.young', 'age.old', 'control.type',
                                   'chron.control.id')
    } else {
      #  If there is no chron-control table, return an empty control.table element.
      control.table <- data.frame(depth = NA,
                                  thickness = NA,
                                  age = NA, 
                                  age.young = NA,
                                  age.old = NA, 
                                  control.type = NA,
                                  chron.control.id = NA)
    }
    
    if ('Default' %in% names(aa[[1]])) {
      meta.table <- data.frame(default     = aa[[1]]$Default,
                               name        = aa[[1]]$ChronologyName,
                               age.type    = aa[[1]]$AgeType,
                               age.model   = aa[[1]]$AgeModel,
                               age.older   = aa[[1]]$AgeOlder,
                               age.younger = aa[[1]]$AgeYounger,
                               chron.id    = aa[[1]]$ChronologyID,
                               date        = aa[[1]]$DatePrepared)
    } else {
      # Return an empty data.frame.
      meta.table <- data.frame(default     = NA,
                               name        = NA,
                               age.type    = NA,
                               age.model   = NA,
                               age.older   = NA,
                               age.younger = NA,
                               chron.id    = NA,
                               date        = NA)
    }
    
    if ('datasets' %in% names(aa[[1]])) {
      parent <- do.call(rbind, aa[[1]]$dataset)
      colnames(parent) <- c('dataset.type', 'dataset.id')
      parent$dataset.id <- as.numeric(parent$dataset.id)
    } else {
      parent <- data.frame(dataset.type = NA,
                           dataset.id = NA)
    }
  }
  
  output <- list(chron.control = control.table,
                 meta = meta.table,
                 parent = parent)
  class(output) <- c('chroncontrol', 'list')
  
  output
            
}
#' @title Function to return chronological control tables from a \code{download} object.
#' @description Using a \code{download}, return the default chron-control table as a \code{data.frame}.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr GET content
#' @param x A single \code{download} object.
#' @param add Should the \code{chroncontrol} be added to the download object (default \code{FALSE})
#' @param verbose logical; should messages on API call be printed?
#' @export
get_chroncontrol.download <- function(x, verbose = TRUE, add = FALSE) {
  chron_id <- x$sample.meta$chronology.id[1]
  
  if (is.na(chron_id)) {
    output <- list(chron.control = data.frame(depth = NA,
                                               thickness = NA,
                                               age = NA, 
                                               age.young = NA,
                                               age.old = NA, 
                                               control.type = NA,
                                               chron.control.id = NA),
                    meta   = data.frame(default     = NA,
                                            name        = NA,
                                            age.type    = NA,
                                            age.model   = NA,
                                            age.older   = NA,
                                            age.younger = NA,
                                            chron.id    = NA,
                                            date        = NA),
                    parent = data.frame(dataset.type = NA,
                                         dataset.id = NA))
    
    class(output) <- c('chroncontrol', 'list')
    
    warning('The download has no assigned chronology id.  Returning an empty data.frame')
    
    return(output)
    
  } else {  
    
    output <- get_chroncontrol(x$sample.meta$chronology.id[1], verbose)
    class(output) <- c('chroncontrol', 'list')
    
    if (names(x$chronologies)[1] == x$sample.meta$chronology.name[1] & add == TRUE) {
      # If we're adding the chronology to the download object:
      x$chronologies[[1]] <- list(chronology = x$chronologies[[1]],
                                  chroncontrol = output)
      return(x)
    }
    
    return(output)
    
  }

}

#' @title Function to return chronological control tables from a \code{download_list} object.
#' @description Using a \code{download_list}, return the default chron-control table as a \code{data.frame}.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr GET content
#' @param x A \code{download_list} object.
#' @param add Should the \code{chroncontrol} be added to the download object (default \code{FALSE})
#' @param verbose logical; should messages on API call be printed?
#' @export
get_chroncontrol.download_list <- function(x, verbose = TRUE, add = FALSE) {
  
  output <- lapply(x, function(y)get_chroncontrol(y, verbose, add = add))
  if (add == FALSE) {
    class(output) <- c('chroncontrol_list', 'list')
  } else {
    class(output) <- c('download_list', 'list')
  }
  return(output)

}
