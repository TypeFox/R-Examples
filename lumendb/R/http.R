#' @title Lumen Database HTTP Requests
#' @description This is the workhorse function for executing API requests for the Lumen Database.
#' @details This is mostly an internal function for executing API requests. In almost all cases, users do not need to access this directly.
#' @param verb A character string containing an HTTP verb, defaulting to \dQuote{GET}.
#' @param path A character string with the API endpoint (should begin with a slash).
#' @param query A list specifying any query string arguments to pass to the API.
#' @param body A character string of request body data.
#' @param base A character string specifying the base URL for the API.
#' @param token A character string containing a Lumen Database API token. If missing, defaults to value stored in environment variable \env{LUMEN_TOKEN}.
#' @param ... Additional arguments passed to an HTTP request function, such as \code{\link[httr]{GET}}.
#' @return A list.
#' @export
lumenHTTP <- function(verb = "GET",
                      path = "", 
                      query = list(),
                      body = "",
                      base = "https://lumendatabase.org",
                      token = Sys.getenv("LUMEN_TOKEN"),
                      ...) {
    url <- paste0(base, path)
    h <- httr::add_headers("Accept" = "application/json", 
                           "AUTHENTICATION_TOKEN" = token)
    if (!length(query)) query <- NULL
    if (verb == "GET") {
        r <- httr::GET(url, query = query, h, ...)
    } else if (verb == "POST") {
        r <- httr::POST(url, body = body, query = query, h, ...)
    } 
    httr::warn_for_status(r)
    return(httr::content(r, "parsed"))
}

#' @title Get a Notice
#' @description Get a Lumen Database Notice by ID
#' @details This function retrieves the details of a Lumen Database notice based upon its Notice ID number. To search for notices, use \code{\link{ldsearch}} instead.
#' @param notice A single numeric value specifying a notice ID.
#' @param ... Additional arguments passed to \code{\link{lumenHTTP}}.
#' @return A list of class \dQuote{lumen_notice}. A \code{summary} method will display some essential details of the notice.
#' @examples
#' \dontrun{
#' n <- ldnotice(1)
#' summary(n)
#' }
#' @export
ldnotice <- function(notice, ...) {
    x <- lumenHTTP(path = paste0("/notices/", notice, ".json"), ...)
    structure(lapply(x, `class<-`, "lumen_notice")[[1]])
}

#' @export
summary.lumen_notice <- function(object, ...) {
    cat(object$type, " notice (", object$id, "): ", object$title, "\n", sep = "")
    cat("Date Received: ", object$date_received, "\n", sep = "")
    cat("Sender:        ", object$sender_name, "\n", sep = "")
    cat("Principal:     ", object$principal_name, "\n", sep = "")
    cat("Recipient:     ", object$recipient_name, "\n", sep = "")
    cat("# of Infringing URLs:  ", sum(sapply(object$works, function(a) length(a$infringing_urls))), "\n", sep = "")
    cat("# of Copyrighted URLs: ", sum(sapply(object$works, function(a) length(a$copyrighted_urls))), "\n\n", sep = "")
    invisible(object)
}

#' @title Get Topics
#' @description Get a list of Lumen Database Topics
#' @details This function retrieves a list of \dQuote{topics} used by the Lumen Database.
#' @param ... Additional arguments passed to \code{\link{lumenHTTP}}.
#' @return A list of objects of class \dQuote{lumen_topic}. The default \code{print} method will display some essential details of each topic.
#' @examples
#' \dontrun{
#' x <- ldtopics()
#' x
#' }
#' @export
ldtopics <- function(...) {
    x <- lumenHTTP(path = paste0("/topics.json"), ...)
    structure(lapply(x[[1]], `class<-`, "lumen_topic"))
}
#' @export
print.lumen_topic <- function(x, ...) {
    cat("Topic (", x$id, "): ", x$name, "\n", sep = "")
    invisible(x)
}


#' @title Search Notices
#' @description Search for Lumen Database Notices
#' @details This function retrieves a list of notices matching a query. Results are paginated by the \code{page} and \code{per_page} arguments. Individual notices can instead be retrieved by their ID using \code{\link{ldnotice}}.
#' @param query A list specifying search query parameters. A reasonable default query would be \code{query = list(term = "joe")} to search for entities containing the word \dQuote{joe}. See \href{https://github.com/berkmancenter/lumendatabase/blob/dev/doc/api_documentation.mkd}{API Documentation} for details of search query terms. 
#' @param page A numeric value specifying which page of results to return. Pagination details are stored in the \code{meta} attribute of the response object.
#' @param per_page A numeric value specifying the number of entities to return in one page. Pagination details are stored in the \code{meta} attribute of the response object.
#' @param verbose A logical (\code{TRUE}, by default) specifying whether to print pagination details to the console.
#' @param ... Additional arguments passed to \code{\link{lumenHTTP}}.
#' @return A list of objects of class \dQuote{lumen_notice}. The default \code{print} method will display some essential details of each topic.
#' @examples
#' \dontrun{
#' # find YouTube-related notices
#' x <- ldsearch(query = list(term = "youtube"))
#' str(x, 1)
#' }
#' @export
ldsearch <- function(query = list(), page = 1, per_page = 10, verbose = TRUE, ...) {
    x <- lumenHTTP(path = paste0("/notices/search"), 
                   query = c(query, list(page = page, per_page = per_page)), ...)
    if (verbose) {
        message(sprintf("Page %s of %s Returned. Response contains %s of %s %s. ", 
                        x$meta$current_page, x$meta$total_pages, 
                        x$meta$per_page, x$meta$total_entries, 
                        ngettext(x$meta$total_entries, "notice", "notices")))
    }
    structure(lapply(x$notices, `class<-`, "lumen_notice"), 
              meta = x$meta,
              class = "lumen_search")
}

#' @export
print.lumen_search <- function(x, ...) {
    for(i in seq_along(x)) {
        cat(paste0("[[", i, "]]\n"))
        summary(x[[i]])
    }
    invisible(x)
}


#' @title Get Entities
#' @description Get a list of Lumen Database entities matching a query
#' @details This function retrieves a list of \dQuote{entities} named in the Lumen Database that match a query. See \href{https://github.com/berkmancenter/lumendatabase/blob/dev/doc/api_documentation.mkd}{API Documentation} for details. Results are paginated by the \code{page} and \code{per_page} arguments.
#' @param query A list specifying search query parameters. A reasonable default query would be \code{query = list(term = "joe")} to search for entities containing the word \dQuote{joe}.
#' @param page A numeric value specifying which page of results to return. Pagination details are stored in the \code{meta} attribute of the response object.
#' @param per_page A numeric value specifying the number of entities to return in one page. Pagination details are stored in the \code{meta} attribute of the response object.
#' @param verbose A logical (\code{TRUE}, by default) specifying whether to print pagination details to the console.
#' @param ... Additional arguments passed to \code{\link{lumenHTTP}}.
#' @return A list of objects of class \dQuote{lumen_entity}. The default \code{print} method will display some essential details of each topic.
#' @examples
#' \dontrun{
#' # return entities matching "joe"
#' ldentities(query = list(term = "joe"))
#' 
#' # use non-default pagination arguments
#' ldentities(query = list(term = "pub"), page = 3, per_page = 5)
#' }
#' @export
ldentities <- function(query = list(), page = 1, per_page = 10, verbose = TRUE, ...) {
    x <- lumenHTTP(path = paste0("/entities/search"), 
                   query = c(query, list(page = page, per_page = per_page)), ...)
    if (verbose) {
        message(sprintf("Page %s of %s Returned. Response contains %s of %s %s. ", 
                        x$meta$current_page, x$meta$total_pages, 
                        ifelse(x$meta$per_page < x$meta$total_entries, 
                               x$meta$per_page, 
                               x$meta$total_entries), 
                        x$meta$total_entries, 
                        ngettext(x$meta$total_entries, "entity", "entities")))
    }
    structure(lapply(x$entities, `class<-`, "lumen_entity"), 
              meta = x$meta,
              class = "lumen_entities")
}

#' @export
print.lumen_entity <- function(x, ...) {
    cat("Entity (", x$id, "): ", x$name, "\n", sep = "")
    invisible(x)
}

# ldcreate <- function(query = list(), ...) {
#    lumenHTTP("POST", path = paste0("/entities/search"), query = query, ...)
# }
