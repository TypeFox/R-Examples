##' Stay authenticated.
##'
##' The auth store keeps the session token after authentication so that all
##' API calls can use it. We can also store other things in it.
##'
##' @format An environment.
##' @keywords internal
session_store <- NULL
makeSessionStore <- function () {
    session_store <<- new.env(hash = TRUE, parent = emptyenv())
    session_store$.globals <- list(prompt=getOption("prompt"))
}
makeSessionStore()

##' Kill the active Crunch session
##' @export
logout <- function () {
    if (is.authenticated()) try(crGET(rootURL("logout")), silent=TRUE)
    deleteSessionInfo()
    options(prompt = session_store$.globals$prompt)
}

##' @importFrom httpcache clearCache
deleteSessionInfo <- function () {
    rm(list=setdiff(ls(envir=session_store), ".globals"), envir=session_store)
    clearCache()
}

##' Authenticate with the Crunch API
##'
##' Note that you can store your Crunch account info in your .Rprofile under
##' "crunch.email" and "crunch.pw" for convenience. If you do so, you can simply
##' \code{login()} to authenticate. For running batch jobs, this could be
##' particularly useful. However, be warned that storing your
##' password in a plain text file such as .Rprofile is a security risk (though
##' perhaps less so than in every .R script you write), and we
##' cannot officially recommend that you do so.
##'
##' If a password is not supplied (or, if no arguments are supplied and only
##' the \code{crunch.email} is specified in .Rprofile), and you are in an
##' interactive session, you will be prompted to enter your password. At
##' present, this is the most secure practice as your password is not stored
##' locally.
##'
##' @param email the email address associated with the user's Crunch account
##' @param password the password associated with the user's Crunch account
##' @param ... additional parameters passed in the authentication. Not
##' currently supported by the Crunch API.
##' @export
login <- function (email=getOption("crunch.email"),
                   password=getOption("crunch.pw"), ...) {
    logout()
    auth <- crunchAuth(email=email, password=password, ...)

    warmSessionCache()

    message("Logged into crunch.io as ", email)
    options(prompt = paste("[crunch]", session_store$.globals$prompt))
    ## Return a Session object
    invisible(session())
}

##' Get various catalogs for your Crunch session
##' @return a list.
##' @export
session <- function () {
    list(datasets=datasetCatalog(),
        projects=ProjectCatalog(crGET(sessionURL("projects"))))
}

crunchAuth <- function (email, password=NULL, ...) {
    ## Validate authentication inputs and then POST to the API
    if (is.null(email)) {
        halt("Must supply the email address associated with your crunch.io account")
    }
    if (is.null(password)) {
        if (interactive()) {
            cat(paste0("Crunch.io password for ", email, ": "))
            unix.alike <- .Platform$OS.type == "unix"
            if (unix.alike) {
                ## Don't print the password being typed
                system("stty -echo")
            }
            password <- readline()
            if (unix.alike) {
                system("stty echo")
                cat("\n")
            }
        } else {
            halt("Must supply a password")
        }
    }

    crPOST(absoluteURL("public/login/", getOption("crunch.api")),
        body=toJSON(list(email=email, password=password, ...)),
        status.handlers=list(`401`=function (response, user=email) {
            halt(paste("Unable to authenticate", user))
        }))
}

jupyterLogin <- function (token) {
    ## Add an auth token as a cookie manually, rather than from a Set-Cookie
    ## response header.
    ## Also modify the user-agent to include "Jupyter"
    set_config(c(config(cookie=paste0("token=", token)),
        add_headers(`user-agent`=crunchUserAgent("jupyter.crunch.io"))))
    warmSessionCache()
}

warmSessionCache <- function () {
    session_store$root <- getAPIroot()
    # session_store$user <- getUser()
    updateDatasetList()
}

is.authenticated <- function () !is.null(session_store$root)
