##############################
## GENRAL UTILITY FUNCTIONS ##
##############################
##Functions in this file:
## stCheckClass       EX:ok
## stCheckFields      EX:ok
## stCheckObs         EX:ok
## stCheckCovars      EX:ok
## stCheckSTcovars    EX:ok
## convertCharToDate  EX:ok
## defaultList        EX:ok

#######################
## Testing functions ##
#######################
##' Test if an object belongs to given class(es).
##'
##' Test if an object belongs to given class(es), and produce reasonable error
##' message if not.
##'
##' @param x Object to test.
##' @param what A character vector naming classes.
##' @param name Character string to be pasted into the error message describing
##'    \code{x}.
##' @return Nothing
##'
##' @examples
##'   ##create a basic object
##'   x <- 1
##'   class(x) <- "test"
##'   ## should be ok
##'   stCheckClass(x, "test", "x")
##'   ## this fails
##'   try( stCheckClass(x, "other", "x") )
##' 
##' @seealso Similar to \code{\link[base:inherits]{inherits}}
##' 
##' @author Johan Lindström
##' @family object checking utilities
##' @export
stCheckClass <- function(x, what, name="Object"){
  if( !inherits(x,what) ){
    stop( paste(name," must belong to one of class(es) ",paste(what,collapse=", ")) )
  }
  return(invisible())
}##function stCheckClass

##' Test if named fields exist in \code{name(x)}, if not the function
##' fails with a suitable error message.
##'
##' @title Test if fields exist in an object.
##' @param x Object to test.
##' @param what A character vector naming that should occur in \code{names(x)}.
##' @param name Character string to be pasted into the error message describing
##'    \code{x}.
##' @return Nothing
##' 
##' @examples
##'   ##load data
##'   data(mesa.model)
##'   ##names present in dta
##'   names(mesa.model$locations)
##' 
##'   ##check for some names
##'   stCheckFields(mesa.model$locations, c("ID","x","lat"))
##'   ##check for non-existant names
##'   try( stCheckFields(mesa.model$locations, c("ID","x","test")) )
##' 
##' @author Johan Lindström
##' @family object checking utilities
##' @export
stCheckFields <- function(x, what, name="Object"){
  if( !all(what %in% names(x)) ){
    stop( paste("Field(s):", paste(what[!(what %in% names(x))],collapse=", "),
                "- missing from", name) )
  }
  return(invisible())
}##function stCheckFields

##' Checks that a observation data.frame is valid.
##'
##' A valid observation data.frame needs to fullfill:
##' \itemize{
##'   \item{Contains fields \code{obs}, \code{date}, and \code{ID}}
##'   \item{All elements in \code{obs$obs} are finte}
##'   \item{\code{obs$date} is one of \code{Date}, \code{numeric},
##'         or \code{integer}}
##'   \item{\code{obs$ID} is \code{character}}
##'   \item{No duplicated observations (same \code{ID} and \code{date})}
##' }
##'
##' @title Check an \code{obs} data.frame.
##' @param obs \code{data.frame} to be checked.
##' @return Nothing
##' 
##' @examples
##'   ##load data
##'   data(mesa.model)
##' 
##'   ##check observations 
##'   stCheckObs( mesa.model$obs )
##'   ##some possible failures
##'   mesa.model$obs <- rbind(mesa.model$obs, mesa.model$obs[1,])
##'   try( stCheckObs( mesa.model$obs ) )
##'   mesa.model$obs$obs[1] <- NaN
##'   try( stCheckObs( mesa.model$obs ) )
##'   mesa.model$obs$date <- as.character( mesa.model$obs$date )
##'   try( stCheckObs( mesa.model$obs ) )
##'   mesa.model$obs$date <- NULL
##'   try( stCheckObs( mesa.model$obs ) )
##' 
##' @author Johan Lindström
##' @family object checking utilities
##' @export
stCheckObs <- function(obs){
  ##should be a data.frame
  stCheckClass(obs, "data.frame", name="obs")
  ##should have certain fields
  stCheckFields(obs, c("obs","date","ID"), name="obs")
  ##check that obs$date is of a suitable type.
  stCheckClass(obs$date, c("integer","numeric","Date"), name="obs$date")
  ##check that obs is of suitable type and finte
  if( !is.numeric(obs$obs) || any(!is.finite(obs$obs)) ){
    stop("Some observation(s) are not finte/numeric.")
  }
  ##check that ID is of suitable type
  stCheckClass(obs$ID, "character", name="obs$ID")
  ##check for duplicates
  if( anyDuplicated(obs[,c("date","ID")]) ){
    stop("Duplicated observations, i.e. multiple elements with same obs$date and obs$ID.")
  }
  
  return(invisible())
}##function stCheckObs
    
##' Checks that data.frame of covariates is valid, making sure that all
##' locations specified in \code{ID.unique} exist. The function will attempt
##' to name each row in \code{covars} using 1) \code{covars$ID}, 2)
##' \code{rownames(covars)}, and 3) \code{as.character(1:dim(covars)[1])}.
##' The field \code{covars$ID} is added if missing and \code{rownames} are removed.
##'
##' @title Check a data.frame of Covariates
##' @param covars data.frame containing covariates, to be checked.
##' @param ID.unique vector with unique IDs that HAVE to be present in the
##'   covariates, typically the observation locations.
##' @return Updated \code{covars} data.frame.
##' 
##' @examples
##'   ##load data
##'   data(mesa.model)
##' 
##'   ##check covariates
##'   tmp <- stCheckCovars( mesa.model$locations, mesa.model$locations$ID )
##'   str(tmp)
##'   ##require non-existant site
##'   try( stCheckCovars( mesa.model$locations, "Bad.Site" ) )
##'   ##drop the ID
##'   mesa.model$locations$ID <- NULL
##'   tmp <- stCheckCovars( mesa.model$locations )
##'   ##ID:s infered from rownames (1-25)
##'   str(tmp)
##' 
##' @author Johan Lindström
##' @family object checking utilities
##' @export
stCheckCovars <- function(covars, ID.unique=character(0)){
  ##should be a data.frame
  stCheckClass(covars, "data.frame", name="covars")

  ##attempt to extract location names from covars$ID or rownames
  if( is.null(covars$ID) ){
    covars$ID <- rownames(covars)
  }
  if( is.null(covars$ID) ){
    warning("No covars$ID or rownames(covars), using ID <- as.character(1:dim(covars)[1])")
    covars$ID <- as.character(1:dim(covars)[1])
  }
  ##make sure that ID:s are characters
  covars$ID <- as.character(covars$ID)
  ##check for duplicates
  if( anyDuplicated(covars$ID) ){
    stop("Duplicated location ID:s in covars$ID / rownames(ID)")
  }
  ##check that we have all the locations we need
  if( !all(ID.unique %in% covars$ID) ){
    stop( paste("Missing locations ",
                paste(ID.unique[!(ID.unique %in% covars$ID)],
                      collapse=", "), "from covar$ID") )
  }
  ##drop rownames
  rownames(covars) <- NULL
  return( covars )
}##function stCheckCovars

##' Checks that array/list of spatio-temporal covariates is valid, making sure
##' that at least all locations specified in \code{ID.unique} exist. The function will
##' attempt to name extract locations ID's from \code{colnames(ST)} and
##' observation dates from \code{rownames(ST)} (using
##' \code{\link{convertCharToDate}}).
##'
##' @title Check an Array/List of Spatio-Temporal Covariates
##' @param ST A 3D-\code{array} containing the ST-covariates, or a
##'   \code{list} of \code{array}:s, the list elements have to be of matching
##'   sizes and have the same \code{rownames} and \code{colnames}; list elemets
##'   are stacked to form a 3D-array.
##' @param date.unique vector with unique dates/times that HAVE to be
##'   present in the ST-covariates, typically the observation time-points.
##' @param ID.unique vector with unique IDs that HAVE to be present in the
##'   ST-covariates, typically the observation locations and un-observation
##'   locations for predictions
##' @return Updated \code{ST} array
##' 
##' @examples
##'   ##load data
##'   data(mesa.model)
##' 
##'   ##check covariates
##'   tmp <- stCheckSTcovars( mesa.model$ST.all, mesa.model$locations$ID )
##'   str(tmp)
##'   ##require non-existant site
##'   try( stCheckSTcovars( mesa.model$ST.all, "Bad.Site" ) )
##'   ##require non-existant site
##'   try( stCheckSTcovars( mesa.model$ST.all, date.unique=1 ) )
##' 
##' @author Johan Lindström
##' @family object checking utilities
##' @export
stCheckSTcovars <- function(ST, ID.unique=character(0), date.unique=integer(0)){
  if( is.null(ST) ){
    ##no ST-covariates
    return(ST)
  }
  ##list, start by collapsing to an array
  if( is.list(ST) ){
    ##extract size, rownames and columnnames
    sz <- dim(ST[[1]])
    col.names <- colnames(ST[[1]])
    row.names <- rownames(ST[[1]])
    ##temporary structure
    tmp <- array(NA, c(sz,length(ST)))
    for(i in 1:length(ST)){
      if( any(dim(ST[[i]])!=sz) || is.null(dim(ST[[i]])) ){
        stop("All elements of list 'ST' must have the same dimensions.")
      }
      if( any(rownames(ST[[i]])!=row.names) || any(colnames(ST[[i]])!=col.names) ){
        stop("row and column names must be equal for all elements of list 'ST'")
      }
      tmp[,,i] <- ST[[i]]
    }#for(i in 1:length(ST))
    dimnames(tmp) <- list(row.names, col.names, names(ST))
    ST <- tmp
  }else{
    ##should be an array
    stCheckClass(ST, "array", name="ST")
    ##and have three dimmensions
    if( length(dim(ST))!=3 ){
      stop("length(dim(ST))!=3")
    }
  }#if( is.list(ST) ){...}else{...}

  ##ID and third dimension names for our spatio-temporal covariate
  if( is.null(colnames(ST)) ){
    warning("No colnames(ST), using as.character(1:dim(ST)[2])")
    colnames(ST) <- as.character(1:dim(ST)[2])
  }
  if( is.null(dimnames(ST)[[3]]) ){
    warning("No dimnames(ST)[[3]], using as.character(1:dim(ST)[3])")
    dimnames(ST)[[3]] <- as.character(1:dim(ST)[3])
  }
  ##rownames should match to Dates, this is a bit harder.
  date <- convertCharToDate( rownames(ST) )
  if( is.null(date) ){
    warning("Unable to coerce rownames(ST) to Date, using date <- 1:dim(ST)[1]")
    date <- 1:dim(ST)[1]
    rownames(ST) <- as.character(date)
  }
  ##check for duplicates
  if( anyDuplicated(colnames(ST)) ){
    stop("Duplicated location, i.e. colnames(ST)")
  }
  if( anyDuplicated(rownames(ST)) ){
    stop("Duplicated dates, i.e. rownames(ST)")
  }
  if( anyDuplicated(dimnames(ST)[[3]]) ){
    stop("Duplicated covariate names, i.e. dimnames(ST)[[3]]")
  }

  ##check that we have all the locations we need
  if( !all(ID.unique %in% colnames(ST)) ){
    stop( paste("Missing location(s)",
                paste(ID.unique[!(ID.unique %in% colnames(ST))],
                      collapse=", "), "from colnames(ST)") )
  }
  ##check that we have all the dates
  if( !all(date.unique %in% date) ){
    stop( paste("Missing date(s)",
                paste(date.unique[!(date.unique %in% date)],
                      collapse=", "), "from rownames(ST)") )
  }
  
  return( ST )
}##function stCheckSTcovars


#######################
## Utility functions ##
#######################
##' Attempts to convert input vector to Date, if that fails tries to convert to
##' double. If conversion induces \code{NA} the function returns \code{NULL}
##' indicating a failure.
##'
##' @title Convert Character to Dates
##' @param x character vector to convert to dates
##' @return a vector of dates, or of doubles or \code{NULL}.
##'
##' @example Rd_examples/Ex_convertCharToDate.R
##' 
##' @author Johan Lindström
##' @family utility functions
##' @export
convertCharToDate <- function(x){
  if( is.null(x) ){
    return(x)
  }
  date <- try( as.Date(x), silent=TRUE)
  if( class(date)!="try-error" && all(!is.na(date)) ){
    ##ok, return value
    return(date)
  }
  suppressWarnings( date <- as.double(x) )
  if( !is.null(date) && all(!is.na(date)) ){
    ##ok, return value
    return(date)
  }
  ##unable to create date vector
  return(NULL)
}

##' Given two lists elemets (by name) missing from the first are copied from the
##' second list (if present). Is used to create default lists, ensuring that all
##' elements expected in the list are present with reasonable values (if not
##' user specified).
##'
##' @title Add Default Elements to Incomplete list
##' @param x A list
##' @param prototype A list with named elements, any elements missing from
##'   \code{x} are replaced with corresponding elements from \code{prototype}.
##' @return Updated version of \code{x}
##'
##' @examples
##'  defaultList(list(a=1,b=4), list(a=3,c="a",d=4))
##' 
##' @author Johan Lindström
##' @family utility functions
##' @export
defaultList <- function(x, prototype=list()){
  if( !is.list(x) ){
    stop("'x' should be a list.")
  }
  if( !is.list(prototype) || is.null(names(prototype)) ){
    stop("'prototype' should be a named list.")
  }
  for(i in names(prototype)){
    if( is.null(x[[i]]) ){
      ##add missing default to x
      x[[i]] <- prototype[[i]]
    }
  }##for(i in names(prototype))
  return(x)
}##function defaultList
