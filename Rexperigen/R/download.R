#' Makes a download request from the server.
#'
#' This function downloads the results of an experiment.
#' 
#' @param sourceURL The source URL for the experiment
#' @param experimentName The experiment name as set in \code{settings.js}
#' @param destination The file to download. By default, all uploaded results are
#' saved in \code{default.csv}, which is the default file to download. 
#' @param auth Whether authentication is needed
#'
#' @return The downloaded data set as data frame
#'
#' @family download
#'
#' @examples
#' \dontrun{
#' downloadExperiment("https...s3.amazonaws.com.myexperiment.index.html", "testing1")
#' downloadExperiment("https...s3.amazonaws.com.myexperiment.index.html", "running", auth = TRUE)
#' }
#' @export
#' 
downloadExperiment <- function(sourceURL, experimentName,
                                destination = "default.csv",
                               auth = FALSE){
    ## if versionMain > 2, we were streaming JSON, so we'll have to convert
    if(versionMain() >= 2){
        request <- checkAuthentication("streamresults", auth)
        res <- API.request(request = request$request,
                   params = list(
                       sourceurl = sourceURL,
                       experimentName = experimentName,
                       file = destination,
                       ndjson = "true"
                   ),
                   auth = request$auth)
        jsonlite::stream_in(textConnection(res), verbose = FALSE)
    }
    else {
        request <- checkAuthentication("makecsv", auth)
        res <- API.request(request = request$request,
                           params = list(
                               sourceurl = sourceURL,
                               experimentName = experimentName,
                               file = destination
                           ),
                           auth = request$auth)
        read.table(text = res, header = TRUE)
    }
}

#' Returns the list of destination files for an experiment.
#'
#' @param sourceURL The source URL for the experiment
#' @param experimentName The experiment name as set in \code{settings.js}
#' @param auth Whether authentication is needed
#'
#' @return The list of destinations
#'
#' @family download
#'
#' @examples
#' \dontrun{
#' getDestinations("https...s3.amazonaws.com.myexperiment.index.html", "running", auth = TRUE)
#' }
#' @export
getDestinations <- function(sourceURL, experimentName, auth = FALSE){
    request <- checkAuthentication("destinations", auth, 2)
    res <- API.request(request = request$request,
                       params = list(
                           sourceurl = sourceURL,
                           experimentName = experimentName
                       ),
                       auth = request$auth)
    tryCatch(
        retval <- jsonlite::fromJSON(res),
        error = function(cond){
            stop(res)
        }
    )
    retval
}

#' Requests the table of users from the server.
#'
#' @param sourceURL The source URL for the experiment
#' @param experimentName The experiment name as set in \code{settings.js}
#' @param auth Whether authentication is needed
#'
#' @return The table of users
#'
#' @family download
#'
#' @examples
#' \dontrun{
#' getUsers("https...s3.amazonaws.com.myexperiment.index.html", "running", auth = TRUE)
#' }
#'
#' @export
getUsers <- function(sourceURL, experimentName, auth = FALSE){
    request <- checkAuthentication("users", auth)
    res <- API.request(request = request$request,
                       params = list(
                           sourceurl = sourceURL,
                           experimentName = experimentName
                       ),
                       auth = request$auth)
    read.table(text = res, header = TRUE)
}

