##' Return base URL for the API calls.
##'
##' Defaults to use beta URL. Not exported.
##' @title base URL
##' @return string for the URL
##' @param dev Should be the beta version of the API be used?
##' @author Francois Michonneau
idig_url <- function(dev=FALSE) {
    if (dev) {
        "http://beta-search.idigbio.org"
    } else {
        "http://search.idigbio.org"
    }
}

##' Return the version number to use for the API calls.
##'
##' The current default is "v2". Not exported.
##' @title API version
##' @param version optional argument giving the version of the API to use
##' @return string for the version to use
##' @author Francois Michonneau
idig_version <- function(version="v2") {
    stopifnot(identical(typeof(version), "character"))
    version
}

##' Parses output of successful query to return a list.
##'
##' Not exported.
##' @title parse successfully returned request
##' @param req the returned request
##' @return a list
##' @author Francois Michonneau
idig_parse <- function(req) {
    txt <- httr::content(req, as="text")
    if (identical(txt, "")) stop("No output to parse", call. = FALSE)
    jsonlite::fromJSON(txt, simplifyVector=FALSE)
}

##' Checks for HTTP error codes and JSON errors.
##'
##' Part 1 of the error checking process. This part handles HTTP error codes and
##' then calls part 2 which handles JSON errors in the responses. Not exported.
##' @title check HTTP code
##' @param req the returned request
##' @return nothing. Stops if HTTP code is >= 400
##' @author Francois Michonneau
idig_check <- function(req) {
  if (req$status_code >= 400) {
    msg <- idig_parse(req)
    stop("HTTP failure: ", req$status_code, "\n", msg$error, "\n",
         msg$context, "\n", msg$name, call. = FALSE)
  }
  idig_check_error(req)
}

##' Checks for error messages that can be returned by the API in JSON.
##'
##' Part 2 of the error checking process. Checks the JSON response for error 
##' messages and stops if any are found. Not exported.
##' @title Check is the request returned an error.
##' @param req the returned request
##' @return nothing. Stops if request contains an error.
##' @author Francois Michonneau
idig_check_error <- function(req) {
  cont <- httr::content(req)
  if (is.list(cont) && exists("error", cont)) {
    stop(paste("Error: ", cont$error, "\n", sep = ""))
  }
}

##' Internal function for GET requests.
##'
##' Generates a GET request and performs the checks on what is returned. Not 
##' exported.
##' @title internal GET request
##' @param path endpoint
##' @param ... additional arguments to be passed to httr::GET
##' @return the request (as a list)
##' @author Francois Michonneau
idig_GET <- function(path, ...) {
    req <- httr::GET(idig_url(), path=paste(idig_version(), path, sep="/"), ...)
    idig_check(req)
    req
}

##' Internal function for POST requests.
##'
##' Generates a POST request and performs the checks on what is returned. Not 
##' exported.
##' @title internal POST request
##' @param path endpoint
##' @param body a list of parameters for the endpoint
##' @param ... additional arguments to be passed to httr::POST
##' @return the request (as a list)
##' @author Francois Michonneau
##' 
idig_POST <- function(path, body, ...) {

    stopifnot(inherits(body, "list"))
    #stopifnot(exists("rq", body))
    
    # Manually encode so we can use auto_unbox=TRUE, see ticket 
    # https://github.com/iDigBio/ridigbio/issues/3
    json <- jsonlite::toJSON(body, auto_unbox=TRUE)
    req <- httr::POST(idig_url(), path=paste(idig_version(), path, sep="/"),
                      body=json, ...)
    idig_check(req)

    req
}

##' Stub function for validating parameters.
##' 
##' Takes list of inputs named by validation rule eg "number":[2, 3] and returns
##' a vector of strings with any validation errors. If the vector is 0 length, 
##' everything is valid. Not exported.
##' @title validate fields
##' @param inputs list of inputs to validate
##' @return boolean
##' @author Matthew Collins
##' 
idig_validate <- function(inputs){
}

##' Format list of field names and indexes.
##' 
##' Some fields returned by the API contain lists or dicts. This function uses
##' a hard coded list of those fields (they are stored in ES with these types
##' so they are known for indexTerms) to generate a list pretty names and 
##' indexes to the returned data OR to drop the fields if formating to a fixed
##' number of columns is not possible. R syntax note:
##'   l[["a"]][["b"]] == l$a$b == l[[c("a", "b")]]
##' Note: indexes assume that the returned JSON is unlisted() first and that 
##' indexTerms and data lists are packed into a list with the keys "indexTerms"
##' and "data". Not exported.
##' @title format field names
##' @param fields list of field names supplied by user
##' @return a list indexed by field name with two values: "data", a list of 
##' field names from data and "indexTerms", a list of field names from 
##' indexTerms
##' @author Matthew Collins
#idig_field_indexes <- function(fields){
#  # looping is old school but keeps the order of fields similar to user input
#  l = list()
#  for (i in fields) {
#    if (i == "flags" ||
#        i == "recordids" ||
#        i == "mediarecords" ||
#        i == "records") {
#      next
#    } else if (i == "geopoint") {
#      l[[paste0(i, ".lat")]] <- c("indexTerms", "geopoint.lat")
#      l[[paste0(i, ".lon")]] <- c("indexTerms", "geopoint.lon")
#    } else if (substr(i, 1, 5) == "data."){
#      l[[i]] <- c("data", substring(i, 6))
#    } else {
#      l[[i]] <- c("indexTerms", i)
#    }
#  }
#  l
#  list("names"=c("uuid", "geopoint.lon", "geopoint.lat"), 
#       "indexTerms"=c("uuid", "geopoint"), 
#       "data"=c())
#}