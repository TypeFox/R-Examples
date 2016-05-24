#' Querying/setting yummlyr option
#'
#' To list all \code{yummlyr} options, just run this function without any parameters provided. To query only one value, pass the first parameter. To set that, use the \code{value} parameter too.
#'
#' The following \code{yummlyr} options are available:
#'
#' \itemize{
#'      \item \code{log}: \code{NULL} or  an optionally passed \emph{logger name} from \pkg{futile.logger} to record all info, trace, debug and error messages.
#'}
#' @param o option name (string). See below.
#' @param value value to assign (optional)
yummlyr_options <- function(o, value) {
    res <- getOption('yummlyr')
    ## just querying
    if (missing(value)) {
        if (missing(o))
            return(res)
        if (o %in% names(res))
            return(res[[o]])
        cat("Possible `yummlyr` options:")
        print(names(res))
        stop("Wrong option queried.")
    } else {
        if (!o %in% names(res))
            stop(paste("Invalid option name:", o))
        ## fix assigning NULL to a list element
        if (is.null(value)) {
            res[o] <- list(NULL)
        } else {
            res[[o]] <- value
        }
        options("yummlyr" = res)
    }
}

#' Process query
#' 
#' Query Yummly API and check return codes
#' @param query string query to execute
perform_query <- function(query) {
    response <- httr::GET(query)
    response_code <- response$status_code
    response_content <- rawToChar(response$content)
    if (response_code == 409) {
        error_massage <- ifelse(grepl("Permission denied", response_content),
                                "Wrong credentials", "API Rate Limit Exceeded")
        stop(error_massage)
    } else if (response_code == 500) {
        stop("Request returned with Internal Server Error, please try again later")
    } else if (response_code == 400) {
        stop("Request is not formatted correctly, please report this error to the developers")
    } else if (response_code != 200) {
        stop(sprintf("Request returned with the following error %s", response_content))
    }
    response_content
}

#' Parse JSONP returned by Yummly for metadata
#' 
#' This function parses JSONP that yummly uses as a response.
#' It is based on assumption that list of elements is returned.
#' @param jsonp jsonp string
parse_jsonp <- function(jsonp) {
    # remove function name and opening parenthesis
    jsonp <- sub('[^\\[|\\{]*', '', jsonp) 
    # remove closing parenthesis
    jsonp <- sub('\\);*$', '', jsonp)
    jsonlite::fromJSON(jsonp)
}

allowed_metadata <- c("allergy", "diet", "ingredient", "cuisine", "course", "holiday")

#' Get metadata
#' 
#' Return information about metadata
#' @param type metadata type
#' @export
get_metadata <- function(type) {
    if (!tolower(type) %in% allowed_metadata) {
        stop(sprintf("Yummly doesn't have any metadata about %s. Allowed metadata: %s",
                     type,
                     paste(allowed_metadata, collapse=", ")))
    }
    metadata[[type]]
}