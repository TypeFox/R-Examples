#' Registers an experiment to the experimenter.
#'
#' Registers the given experiment to the experimenter.
#' Returns the server's response. Do provide the credentials
#' before calling this function
#'
#' @param sourceURL The source URL for the experiment
#' @param experimentName The experiment name as set in \code{settings.js}
#'
#' @return The server response.
#'
#' @family registration
#'
#' @examples
#' \dontrun{
#' registerExperiment("https...s3.amazonaws.com.myexperiment.index.html", "running")
#' }
#' 
#' @export
registerExperiment <- function(sourceURL, experimentName){
    experimenter <- checkLogin()
    API.request(request = "digest/registration",
                params = list(
                    experimenter = experimenter,
                    sourceurl = sourceURL,
                    experimentName = experimentName
                ),
                auth = TRUE,
                method = "POST"
                )
}

#' Removes the registration of the experiment
#'
#' Removes the registration of the experiment to the logged in experimenter.
#' Does not remove any data though! All of the data will be available to
#' anyone.
#'
#' @param sourceURL The source URL for the experiment
#' @param experimentName The experiment name as set in \code{settings.js}
#'
#' @return The server response.
#'
#' @examples
#' \dontrun{
#' removeRegistration("https...s3.amazonaws.com.myexperiment.index.html", "running")
#' }
#' 
#' @family registration
#'
#' @export
removeRegistration <- function(sourceURL, experimentName){
    experimenter <- checkLogin()
    API.request(request = "digest/registration",
                params = list(
                    experimenter = experimenter,
                    sourceurl = sourceURL,
                    experimentName = experimentName
                ),
                auth = TRUE,
                method = "DELETE"
                )
}


#' Get a list of registered experiments.
#'
#' Returns a list of the registered experiments for the logged in
#' experimenter.
#'
#' @return The parsed list of experiments
#'
#' @family registration
#'
#' @examples
#' \dontrun{
#' getRegisteredExperiments()
#' }
#' 
#' @export
getRegisteredExperiments <- function(){
    experimenter <- checkLogin()
    res <- API.request(request = "digest/registration",
                       params = list(experimenter = experimenter),
                       auth = TRUE,
                       method = "GET")
    jsonlite::fromJSON(res)
}


#' Checks whether the user is logged in.
#'
#' Checks whether the user is logged in. It also throws an exception if the
#' server is of an old version and logging in is not supported.
#'
#' @return The experimenter login user name
#'
#' @examples
#' options(Rexperigen.server.version="2.0.0")
#' setExperigenCredentials("my", "credentials", check = FALSE)
#' checkLogin()
#' 
#' @export
checkLogin <- function(){
    experimenter <- getOption("Rexperigen.experimenter")
    if(experimenter == ""){
        stop("Not logged in.")
    }
    if(versionMain() < 2){
        stop("Registering experiments not supported by server.")
    }
    experimenter
}
