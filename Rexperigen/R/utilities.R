#' Server URL preparation.
#'
#' Prepares the server URL by checking for the presence of a protocol
#' prefix, and prefixing \code{http://} if there is none. Also, returns
#' the URL with a trailing \code{/} if that was not present either.
#'
#' @param server A server name.
#' @return A clean URL with protocol prefix and trailing slash.
#' @family development
#' @examples
#' prepare.server.URL("db.phonologist.org")
#' prepare.server.URL("https://myserver.com/")
#'
#' @export
prepare.server.URL <- function(server){
    url <- ifelse(grepl("://", server), server, paste0("http://",server))
    url <- ifelse(substr(url, nchar(url), nchar(url)) == "/", url, paste0(url, "/"))
    url
}

#' Create an URL for a request to the Experigen API.
#'
#' @param server The server name (no need to be clean)
#' @param request The request verb to the server
#' @param params A list of name=value pairs that will be passed
#' to the server
#' @return A URL that can be queried.
#'
#' @family development
#' @examples
#' create.API.request.URL("db.phonologist.org", "users", list(sourceurl="x.y", experimentName="z"))
#' @export
create.API.request.URL <- function(server, request, params = list()){
   url <- paste0(prepare.server.URL(server), request)
    if(length(params) > 0){
        parampairs <- c()
        for(n in names(params)){
            parampairs <- c(parampairs, paste(n, params[n], sep = "="))
        }
        url <- paste(url, paste(parampairs, collapse = "&"), sep = "?")
    }
    url
}

#' Run a request to the API
#'
#' @param server The server name (no need to be clean)
#' @param request The request verb to the server
#' @param params A list of \code{name=value} pairs that will be passed
#' to the server
#' @param auth Do we send authentication? If this is \code{FALSE}, we don't.
#' @param method The HTTP method to use (\code{GET | POST | PUT | DELETE}).
#' 
#' @family development
#'
#' @export
API.request <- function(server = getOption("Rexperigen.server"),
                        request = "version",
                        params = list(),
                        auth = FALSE,
                        method = "GET"){
    url <- create.API.request.URL(server, request, params)
    if(auth){
        res <- tryCatch(
            RCurl::getURL(url,
                          username = getOption("Rexperigen.experimenter"),
                          password = getOption("Rexperigen.password"),
                          httpauth = RCurl::AUTH_DIGEST,
                          customrequest = method),
            error = function(e){
                stop(paste("Error downloading results:",e$message))
            }
        )
        return(res)
    }
    else {
        res <- tryCatch(
            RCurl::getURL(url,
                          customrequest = method),
            error = function(e){
                stop(paste("Error downloading results:",e$message))
            }
        )
        return(res)
    }
}


#' Returns the version of an Experigen server
#'
#' @param server The server name
#' @return A string with the Experigen server version
#' @examples
#' server.version("db.phonologist.org")
#' @family development
#' @export
server.version <- function(server = getOption("Rexperigen.server")){
    url <- prepare.server.URL(server)
    if(RCurl::url.exists(paste0(url,"version"))){
        RCurl::getURL(paste0(url,"version"))
    }
    else {
        if(RCurl::url.exists(paste0(url,"makecsv.cgi"))){
            "1.0.0"
        }
        else {
            NO_SERVER_ERROR
        }
    }
}

#' Error string returned when no server is found
#' @export
NO_SERVER_ERROR <- "no Experigen server there"


#' Checks whether the authentication is supported by the server
#' and handles the different POST request routes.
#'
#' @param request The needed request to be made (without \code{.cgi}).
#' @param auth Whether authentication is requested in R.
#' @param version.needed The minimum version of the request. If it is
#' \code{1}, and the server is version 1, the function will append
#' \code{.cgi} to the request.
#'
#' @examples
#' checkAuthentication("makecsv", TRUE)
#' checkAuthentication("registration", TRUE, 2)
#'
#' @family development
#'
#' @return A list: \enumerate{
#' \item \code{request}: the request to send
#' \item \code{auth}: whether to still send an authenticated request.
#' }
#'
#' @export
checkAuthentication <- function(request, auth, version.needed = 1){
    if(version.needed < 2){
        request <- ifelse(versionMain() < 2, paste0(request, ".cgi"), request)
    }
    if(auth && versionMain() < 2){
        warning("Experiment registration is not supported by the server, requesting without authentication.")
        auth <- FALSE
    }
    if(auth){
        request <- paste0("digest/", request)
    }
    list(request = request,
         auth = auth)
}

#' Returns the cleaned source URL.
#'
#' Experigen has been using some URL cleaning since version 1. This function
#' asks the server (if version >= 2) to return how it would clean the given
#' URL. For older versions, you still have to replace /'s and stuff by .'s
#' by yourself, sorry!
#'
#' @param sourceURL The source URL to clean
#' @return The cleaned URL
#'
#' @examples
#' \dontrun{
#' cleanURL("https://s3.amazonaws.com/myexperiment/index.html")
#' }
#' 
#'
#' @export
cleanURL <- function(sourceURL){
    if(versionMain() < 2){
        stop("This request is not supported by the server.")
    }
    res <- API.request(request = "cleanURL",
                       params = list(sourceurl = URLencode(sourceURL, TRUE)))
    if(length(res) > length(sourceURL)){ # this has got to be an error
        stop("There was a problem with the conversion -- maybe you are using illegal characters? Be on the lookout for -'s and _'s!")
    }
    res
}
