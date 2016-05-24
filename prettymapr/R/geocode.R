#' Geocode Locations
#'
#' Geocode locations using the
#' \href{https://developers.google.com/maps/documentation/geocoding/intro}{Google Web API} or
#' the \href{https://pickpoint.io/}{PickPoint.io API}. Implemented
#' from the \code{ggmap:geocode} function from the \code{ggmap}
#' package (\url{https://cran.r-project.org/package=ggmap})
#' by David Kahle
#'
#' @param location A character vector of locations to pass to the geocoding API.
#' @param output One of \code{data.frame} or \code{list}. If \code{data.frame},
#'   the results are distilled into columns: query, source, status, rank, address, lon, lat, bbox_n,
#'   bbox_e, bbox_s, bbox_w, and id. If \code{list}, the raw JSON output from the geocoding API
#'   is returned as a \code{list} (containing lists). The output of a failed geocode return will
#'   always have a \code{$status} attribute describing the failure.
#' @param source One of "default", "google" or "pickpoint". If "default", the function
#'   calls \code{options("prettymapr.geosource")} or chooses "pickpoint" if none is set.
#'   If using "pickpoint", please \href{https://pickpoint.io/users/sign_up}{sign up for your own (free) API key}
#'   to avoid using the default excessively.
#' @param messaging \code{TRUE} if verbose messaging is desired.
#' @param limit The number of results to return per query. This refers to individual locations, for which
#'   ambiguous queries may return multiple results (e.g. Halifax, Nova Scotia; Halifax, United Kingdom, etc.).
#'   The default is 1. Pass 0 if no limit on queries is desired.
#' @param key API key if \code{source="pickpoint"}.
#' @param sensor \code{TRUE} if the location is generated from a sensor.
#' @param ... A number of key/value pairs to append to the URL, specifying further options specific to each
#'   API. Google users may wish to provide \code{client} and \code{signature} arguments for use with
#'   the enterprise version with the API, or specify additional constraints on geocoding.
#'
#' @return A \code{list} or \code{data.frame}; see documentation for \code{output} argument.
#'
#' @export
#'
#' @examples
#' #don't test to speed up checking time
#' \donttest{
#' geocode("wolfville, ns")
#' geocode("wolfville, ns", output="list")
#' geocode("halifax", limit=0)
#' geocode("Paddy's Pub Wolfville NS", source="google")
#' geocode(c("Houston, TX", "San Antonio TX", "Cleavland OH"), source="google")
#'
#' #fails quietly
#' geocode("don't even think about geocoding this")
#' geocode("don't even think about geocoding this", output="list")
#' }
#'
geocode <- function(location, output=c("data.frame", "list"), source = "default", messaging = FALSE, limit=1, key=NULL,
                    sensor=FALSE, ...) {

  #source can be "google", "dsk", or "pickpoint"...but dsk doesn't appear to be working
  #output ban be "list" or "data.frame"
  output <- match.arg(output)
  if(!is.logical(messaging)) stop("messaging must be TRUE or FALSE")
  if(!is.logical(sensor)) stop("sensor must be TRUE or FALSE")
  urlargs <- list(...)

  # find source based on input
  if(source=="default") {
    opts <- options("prettymapr.geosource")
    if(!is.null(opts$prettymapr.geosource)) {
      source <- opts$prettymapr.geosource
    }

    if(!(source %in% c("google", "dsk", "pickpoint"))) {
      if(messaging) message("No default source set, using pickpoint")
      source <- "pickpoint"
    }
  }

  #"import" foreach
  foreach <- foreach::foreach
  `%do%` <- foreach::`%do%`

  #deal with vectors input here
  if(length(location) > 10 && is.null(key) && source=="pickpoint") {
    stop("Please get your own PickPoint.io API key at https://pickpoint.io/users/sign_up")
  }

  if(length(location) > 1) {
    loc<-NULL;rm(loc) #trick CMD check
    if(output=="data.frame") {
      return(foreach(loc=location, .combine=rbind) %do% {
        geocode(loc, output, source, messaging, limit, key, sensor, ...)
        })
    } else {
      return(foreach(loc=location) %do% {
        geocode(loc, output, source, messaging, limit, key, sensor, ...)
      })
    }
  }

  # return NA for location == "", NA, or NULL
  if(location == "" || is.na(location) || is.null(location)) {
    return(failedGeocodeReturn(output, source, location))
  }

  # format the url
  if(source=="google" || source=="dsk") {
    urlargs$sensor <- tolower(as.character(sensor))
  }
  if(length(urlargs) > 0) {
    combineargs <- function(...) {paste(..., sep="&")}
    key<-NULL;rm(key) #CMD Check trick
    moreargs <- foreach(key=names(urlargs), .combine=combineargs, .multicombine=TRUE) %do% {
      return(sprintf("%s=%s", key, urlargs[[key]]))
    }
  } else {
    moreargs <- ""
  }

  location4url <- chartr(" ", "+", location)
  posturl <- paste(location, moreargs, sep = "&")

  if(source == "google"){
    url_string <- paste0("http://maps.googleapis.com/maps/api/geocode/json?address=", posturl)
  } else if(source == "dsk"){
    url_string <- paste0("http://www.datasciencetoolkit.org/maps/api/geocode/json?address=", posturl)
  } else if(source=="pickpoint") {
    if(is.null(key)) {
      key <- "yxsN6E8EYYHFF_xsa_uL"
      message("Using default API key, if batch geocoding please get your own (free) API key at https://pickpoint.io/users/sign_up")
    }
    url_string <- paste0("https://api.pickpoint.io/v1/forward?key=", key , "&q=", posturl)
  } else {
    stop("Unknown source provided: ", source)
  }

  url_string <- utils::URLencode(url_string)
  url_hash   <- digest::digest(url_string)

  # lookup info if on file
  if(isGeocodedInformationOnFile(url_hash)){

    if(messaging) message("Using stored information.")
    gc <- get(".GeocodedInformation", envir = .GlobalEnv)[[url_hash]]

  } else {

    if(messaging) message(paste("contacting ", url_string, "...", sep = ""), appendLF = F)
    # message user
    message(paste0("Information from URL : ", url_string))
    message("Queries in this session for ", source, ": ", getGeocodeCount(source), " (limit of 2500/day for both google and pickpoint)")
    # geocode
    connect <- url(url_string)
    lines <- try(readLines(connect, warn = FALSE), silent = TRUE)
    close(connect)

    if(class(lines) == "try-error"){
      warning(
        "  geocoding failed for \"", location, "\".\n",
        "  if accompanied by 500 Internal Server Error with using dsk, try google."
      )
      return(failedGeocodeReturn(output, source, location))
    } else if((length(lines)==1) && (lines=="Unauthorized")) {
      warning("Invalid API specified: ", key)
      return(failedGeocodeReturn(output, source, location))
    }

    gc <- rjson::fromJSON(paste(lines, collapse = ""))
    if(messaging) message(" done.")

    # temporarily save it
    storeGeocodedInformation(url_hash, source, gc)

  }

  #standardize output to query, source, address, lon, lat, bbox_n, bbox_e, bbox_s, bbox_w, id
  #TODO add class, type fields $results[[1]]$types[2], $results[[1]]$types[1] and [[1]]$class, [[1]]$type
  if(length(gc)==0) {
    return(emptyGeocodeReturn(output, source, location))
  }

  if(output=="data.frame") {
    i<-NULL;rm(i) #CMD Check trick
    if(source=="google" || source=="dsk") {
      if(limit==0) {
        limit <- length(gc$results)
      }
      foreach(i=1:min(length(gc$results), limit), .combine=rbind) %do% {
        data.frame(query=location,
                  source=source,
                  status=gc$status,
                  rank=i,
                  address=gc$results[[1]]$formatted_address,
                  lon=as.numeric(gc$results[[i]]$geometry$location$lng),
                  lat=as.numeric(gc$results[[i]]$geometry$location$lat),
                  bbox_n=as.numeric(gc$results[[i]]$geometry$viewport$northeast$lat),
                  bbox_e=as.numeric(gc$results[[i]]$geometry$viewport$northeast$lng),
                  bbox_s=as.numeric(gc$results[[i]]$geometry$viewport$southwest$lat),
                  bbox_w=as.numeric(gc$results[[i]]$geometry$viewport$southwest$lng),
                  id=gc$results[[i]]$place_id, stringsAsFactors = FALSE)
      }
    } else if(source=="pickpoint") {
      if(limit==0) {
        limit <- length(gc)
      }
      foreach(i=1:min(length(gc), limit), .combine=rbind) %do% {
        data.frame(query=location,
                   source=source,
                   status="OK",
                   rank=i,
                   address=gc[[i]]$display_name,
                   lon=as.numeric(gc[[i]]$lon),
                   lat=as.numeric(gc[[i]]$lat),
                   bbox_n=as.numeric(gc[[i]]$boundingbox[2]),
                   bbox_e=as.numeric(gc[[i]]$boundingbox[4]),
                   bbox_s=as.numeric(gc[[i]]$boundingbox[1]),
                   bbox_w=as.numeric(gc[[i]]$boundingbox[3]),
                   id=gc[[i]]$place_id, stringsAsFactors = FALSE)
      }
    }
  } else {
    return(gc)
  }

}

geoInfoDoesntExist <- function(){
  !(".GeocodedInformation" %in% ls(envir = .GlobalEnv, all.names =  TRUE))
}

storeGeocodedInformation <- function(url_hash, source, data){
  .GeocodedInformation <- NULL; rm(.GeocodedInformation) #trick CMD check

  if(geoInfoDoesntExist()) .GeocodedInformation <<- list()

  db <- get(".GeocodedInformation", envir = .GlobalEnv)
  if(is.null(db[[source]])) {
    db[[source]] <- 0
  } else {
    db[[source]] <- db[[source]] + 1
  }

  placesOnFile <- names(db)
  db <- c(db, list(data))
  names(db) <- c(placesOnFile, url_hash)

  .GeocodedInformation <<- db

  invisible()
}

getGeocodeCount <- function(source) {
  if(geoInfoDoesntExist()) return(0)
  db <- get(".GeocodedInformation", envir = .GlobalEnv)
  if(is.null(db[[source]])) {
    0
  } else {
    db[[source]]
  }
}

retrieveGeocodedInformation <- function(url_hash){
  if(geoInfoDoesntExist()) return(NA)
  get(".GeocodedInformation", envir = .GlobalEnv)[[url_hash]]
}

isGeocodedInformationOnFile <- function(url_hash){
  if(geoInfoDoesntExist()) return(FALSE)
  if(!(url_hash %in% names(get(".GeocodedInformation", envir = .GlobalEnv)))) return(FALSE)
  TRUE
}

clearGeocodedInformation <- function(){
  # suppress in case it doesn't exist
  suppressWarnings(rm(".GeocodedInformation", envir = .GlobalEnv))
  invisible()
}

failedGeocodeReturn <- function(output, source, location, status="Failed"){
  out <- list(query=location,
             source=source,
             status=status,
             rank=NA,
             address=NA,
             lon=NA,
             lat=NA,
             bbox_n=NA,
             bbox_e=NA,
             bbox_s=NA,
             bbox_w=NA,
             id=NA)
  if(output == "data.frame") {
    data.frame(out, stringsAsFactors = FALSE)
  } else {
    out
  }
}

emptyGeocodeReturn <- function(output, source, location) {
  failedGeocodeReturn(output, source, location, status="Emtpy")
}


