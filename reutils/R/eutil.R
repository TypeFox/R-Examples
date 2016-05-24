#' @include eutil-error.R
#' @importFrom RCurl basicHeaderGatherer basicTextGatherer curlEscape postForm
#'   getURLContent
NULL

#' Class \code{"eutil"}: Reference classes that hold the response from EUtils
#' requests.
#' 
#' The reference classes \code{\linkS4class{eutil}}, \code{\linkS4class{einfo}},
#' \code{\linkS4class{esearch}}, \code{\linkS4class{esummary}},
#' \code{\linkS4class{efetch}}, \code{\linkS4class{elink}}, \code{\linkS4class{epost}},
#' \code{\linkS4class{egquery}}, \code{\linkS4class{espell}}, and
#' \code{\linkS4class{ecitmatch}} implement the request generator for interaction
#' with the NCBI services.
#' They should not be used direcly, but initialized through the respective
#' constructor functions \code{\link{einfo}}, \code{\link{esearch}},
#' \code{\link{esummary}}, \code{\link{efetch}}, \code{\link{elink}},
#' \code{\link{epost}}, \code{\link{egquery}}, \code{\link{espell}}, and
#' \code{\link{ecitmatch}}.
#' 
#' @field params A named \code{list} of query parameters.
#' @field errors A \code{\linkS4class{eutil_error}} object.
#' @field content Result of an Entrez request stored as a character vector.
#' 
#' @section Extends: All reference classes extend and inherit methods from
#'     \code{"\linkS4class{envRefClass}"}. Furthermore, \code{"einfo"},
#'     \code{"esearch"}, \code{"esummary"}, \code{"efetch"}, \code{"elink"},
#'     \code{"epost"}, \code{"egquery"}, \code{"espell"}, and \code{"ecitmatch"}
#'     all extend the \code{"eutil"} class.
#' 
#' @seealso \code{\linkS4class{eutil}}, \code{\link{einfo}},
#' \code{\link{esearch}}, \code{\link{esummary}}, \code{\link{efetch}},
#' \code{\link{elink}}, \code{\link{epost}}, \code{\link{egquery}},
#' \code{\link{espell}}, and \code{\link{ecitmatch}}.
#' 
#' @name eutil-class
#' @aliases eutil-class einfo-class esearch-class esummary-class efetch-class
#'          elink-class epost-class egquery-class espell-class ecitmatch-class
#' @keywords classes internal
#' @export
#' @examples
#' showClass("eutil")
eutil <- setRefClass(
    Class   = "eutil",
    fields  = list(params = "list", errors = "eutil_error", content = "character"),
    methods = list(
      initialize = function() {
        .self$params  <- list()
        .self$errors  <- eutil_error()
        .self$content <- NA_character_
      },
      ##
      ## public methods
      ##
      xmlValue = function(xpath, as = "character", default = NA_character_) {
        'Extract the text value of XML leaf nodes given a valid XPath expression.'
        xvalue(get_content("xml"), xpath, as, default)
      },
      xmlName = function(xpath, as = "character", default = NA_character_) {
        'Extract the tag names of XML nodes given a valid XPath expression.'
        xname(get_content("xml"), xpath, as, default)
      },
      xmlAttr = function(xpath, name, as = "character", default = NA_character_) {
        'Extract the value of XML attributes given a valid XPath expression
        and an attribute name'
        xattr(get_content("xml"), xpath, name, as, default)
      },
      xmlSet = function(xpath, ...) {
        'Extract a set of XML nodes given a valid XPath expression.'
        xset(get_content("xml"), xpath, ...)
      },
      ##
      ## internal but documented methods
      ##
      get_url = function() {
        "Return the URL used for an Entrez query; should not be used directly,
        use \\code{\\link{getUrl}} instead."
        return(query_url('GET'))
      },
      get_error = function() {
        "Return \\code{\\linkS4class{eutil_error}}s; should not be used directly,
        use \\code{\\link{getError}} instead."
        return(.self$errors)
      },
      get_content = function(as = "text", ...) {
        "Return the results of an Entrez query as text, xml, json, a textConnection,
        or parsed to a list or data.frame; should not be used directly, use
        \\code{\\link{content}} instead."
        as <- match.arg(as, c("text", "xml", "json", "textConnection", "parsed"))
        check_retmode(as)
        switch(as,
          text   = .self$content,
          xml    = savely_parse_xml(.self$content, ...),
          json   = savely_parse_json(.self$content, ...),
          parsed = parse_content(.self, ...),
          textConnection = textConnection(.self$content, ...)
        )
      },
      perform_query = function(method = "GET", ...) {
        "Perform an Entrez query using either http GET or POST requests;
        should not be used directly."
        verbose <- isTRUE(getOption("reutils.verbose.queries"))
        if (verbose) {
          cat("Perfoming an", sQuote(eutil()), "query ...\n")
        }
        method <- match.arg(method, c("GET", "POST"))
        ## update an object with new query parameters
        .email <- getOption("reutils.email")
        if (is.null(.email) || grepl("^Your\\.name\\.here.+", .email, ignore.case = TRUE)) {
          warning("NCBI requests that you provide an email address with each query to their API.\n",
                  "Please set the global option ", sQuote("reutils.email"), " to your address to make",
                  " this message go away.", call. = FALSE, immediate. = FALSE)
        }
        .params <- list(...)
        .params <- compact(Reduce(merge_list, list(.params, params, list(email = .email, tool = "reutils"))))
        .self$params <- .params
        
        opts <- list(connecttimeout = 10)
        hg <- RCurl::basicHeaderGatherer()
        opts$headerfunction <- hg$update
        tg <- RCurl::basicTextGatherer()
        opts$writefunction <- tg$update
        
        if (method == "POST") {
          e <- tryCatch({
            RCurl::postForm(query_url("POST"), .params = .self$params,
                            .opts = opts, style = "POST")
          }, error = function(e) e$message)
        } else if (method == "GET") {
          if (verbose) {
            cat(ellipsize(query_url("GET")), "\n")
          }
          e <- tryCatch({
            RCurl::getURLContent(query_url("GET"), .opts = opts)
          }, error = function(e) e$message)
        }
        .self$content <- as.character(tg$value())
        if (is.null(e) || !nzchar(e)) {
          header <- as.list(tryCatch(hg$value()), error = function(e) {
            c(status = 500, statusMessage = e$message)
          })
          status <- as.numeric(header$status)
          statusmsg <- header$statusMessage
          if (status != 200) {
            .self$errors$error <- paste0("HTTP error: Status ", status, "; ", statusmsg)
            warning(errors$error, call. = FALSE, immediate. = TRUE)
          }
        } else {
          .self$errors$error <- paste0("CurlError: ", e)
          warning(errors$error, call. = FALSE, immediate. = TRUE)
        }
      },
      ##
      ## helper methods (undocumented)
      ##
      query_url = function(method) {
        host <- switch(eutil(),
                       egquery = "http://eutils.ncbi.nlm.nih.gov/gquery",
                       ecitmatch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/ecitmatch.cgi",
                       paste0('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/', eutil(),'.fcgi')
        )
        if (method == "GET") {
          fields <- paste(RCurl::curlEscape(names(.self$params)), RCurl::curlEscape(.self$params),
                          sep = "=", collapse = "&")
          paste0(host, "?", fields)
        } else {
          host
        }
      },
      eutil = function() {
        class(.self)
      },
      database = function() {
        .self$params$db
      },
      rettype = function() {
        .self$params$rettype
      },
      retmode = function() {
        .self$params$retmode
      },
      no_errors = function() {
        .self$errors$all_empty()
      },
      check_retmode = function(as) {
        if (!is.null(retmode())) {
          if (retmode() == 'xml' && as %ni% c("text", "xml", "parsed"))
            stop("Cannot return data of retmode ", dQuote(retmode()), " as ", dQuote(as), ".", call. = FALSE)
          if (retmode() == 'json' && as %ni% c("text", "json", "parsed"))
            stop("Cannot return data of retmode ", dQuote(retmode()), " as ", dQuote(as), ".", call. = FALSE)
          if (retmode() == 'text' && as %ni% c("text", "textConnection"))
            stop("Cannot return data of retmode ", dQuote(retmode()), " as ", dQuote(as), ".", call. = FALSE)
          if (retmode() == 'asn.1' && as %ni% c("text", "textConnection"))
            stop("Cannot return data of retmode ", dQuote(retmode()), " as ", dQuote(as), ".", call. = FALSE)
        }
      }
    )
  )

savely_parse_xml <- function(x, ...) {
  tryCatch(XML::xmlParse(x, asText = TRUE, error = NULL, ...),
           "XMLError" = function(e) {
             errmsg <- paste("XML parse error:", e$message)
             XML::xmlParseString(paste0("<ERROR>", errmsg, "</ERROR>"))
           },
           "error" = function(e) {
             errmsg <- paste("Simple error:", e$message)
             XML::xmlParseString(paste0("<ERROR>", errmsg, "</ERROR>"))
           })
}

savely_parse_json <- function(x, ...) {
  intent <- list(...)$intent %||% 2
  jsonlite::prettify(x, indent = intent)
}

parse_content <- function(object, ...) {
  switch(object$eutil(),
    einfo     = parse_einfo(object, ...),
    esearch   = parse_esearch(object, ...),
    epost     = parse_epost(object, ...),
    esummary  = parse_esummary(object, ...),
    elink     = parse_linkset(object, ...),
    ecitmatch = parse_ecitmatch(object, ...),
    "Not yet implemented"
  )
}

#' Extract the data content from an Entrez request
#' 
#' There are five ways to access data returned by an Entrez request: as a character
#' string \code{(as = "text")}, as a \code{\link{textConnection}}
#' \code{(as = "textConnection")}, as an \code{\linkS4class{XMLInternalDocument}}
#' \code{(as = "xml")} or \code{json} object \code{(as = "json")}
#' (depending on the \code{retmode} with which the request was performed),
#' or parsed into a native R object, e.g. a \code{list} or a \code{data.frame}
#' \code{(as = "parsed")}.
#' 
#' @param x An \code{\linkS4class{eutil}} object.
#' @param as Type of output: \code{"text"}, \code{"xml"}, \code{"json"},
#' \code{"textConnection"}, or \code{"parsed"}. \code{content} attempts to
#' figure out the most appropriate output type, based on the \code{retmode} of
#' the object.
#' @param ... Further arguments passed on to methods.
#' @seealso
#'    \code{\link{einfo}}, \code{\link{esearch}}, \code{\link{esummary}},
#'    \code{\link{efetch}}, \code{\link{elink}}, \code{\link{epost}},
#'    \code{\link{egquery}}, \code{\link{espell}}, \code{\link{ecitmatch}}.
#' @export
#' @examples
#' \dontrun{
#' ## einfo() defaults to retmode 'xml'
#' e <- einfo()
#' 
#' ## automatically return data as an 'XMLInternalDocument'.
#' if (e$no_errors()) {
#'   content(e)
#' 
#'   ## return the XML data as character string.
#'   cat(content(e, "text"))
#' 
#'   ## return DbNames parsed into a character vector.
#'   content(e, "parsed")
#' }
#' 
#' ## return data as a JSON object
#' e2 <- einfo(db = "gene", retmode = "json")
#' if (e2$no_errors()) {
#'   content(e2)
#' }
#' 
#' ## return a textConnection to allow linewise reading of the data.
#' x <- efetch("CP000828", "nuccore", rettype = "gbwithparts", retmode = "text")
#' con <- content(x, as = "textConnection")
#' readLines(con, 2)
#' close(con)
#' }
setGeneric("content", function(x, ...) standardGeneric("content"))
#' @describeIn content
setMethod("content", "eutil", function(x, ...) {
  as <- list(...)$as %||% match.arg(x$retmode(), c("text", "xml", "json"))
  x$get_content(as)
})

#' getError
#' 
#' Retrieve a http or XML parsing error from an \code{\linkS4class{eutil}} object.
#' 
#' @param x An \code{\linkS4class{eutil}} object.
#' @param ... Further arguments passed on to methods.
#' @return An \code{\linkS4class{eutil_error}} object.
#' @seealso
#'    \code{\link{einfo}}, \code{\link{esearch}}, \code{\link{esummary}},
#'    \code{\link{efetch}}, \code{\link{elink}}, \code{\link{epost}},
#'    \code{\link{egquery}}, \code{\link{espell}}, \code{\link{ecitmatch}}.
#' @export
#' @examples
#' \dontrun{
#' e <- efetch("Nonsensical_accession_nr", "protein", rettype = "fasta")
#' getError(e)
#' }
setGeneric("getError", function(x, ...) standardGeneric("getError"))
#' @describeIn getError
setMethod("getError", "eutil", function(x, ...) {
  x$get_error()
})

#' getUrl
#' 
#' Retrieve the URL used to perform an Entrez E-Utilities query.
#' 
#' @param x An \code{\linkS4class{eutil}} object.
#' @param ... Further arguments passed on to methods.
#' @return A character string.
#' @seealso
#'    \code{\link{einfo}}, \code{\link{esearch}}, \code{\link{esummary}},
#'    \code{\link{efetch}}, \code{\link{elink}}, \code{\link{epost}},
#'    \code{\link{egquery}}, \code{\link{espell}}, \code{\link{ecitmatch}}.
#' @export
#' @examples
#' \dontrun{
#' e <- efetch("AV333213.1", "protein", rettype = "fasta")
#' getUrl(e)
#' }
setGeneric("getUrl", function(x, ...) standardGeneric("getUrl"))
#' @describeIn getUrl
setMethod("getUrl", "eutil", function(x, ...) {
  x$get_url()
})


#' performQuery
#' 
#' @param x An \code{\linkS4class{eutil}} object.
#' @param method One of \dQuote{GET} or \dQuote{POST}.
#' @param ... Further arguments passed on to methods.
#' @export
#' @keywords internal
setGeneric("performQuery", function(x, method = "GET", ...) standardGeneric("performQuery"))
#' @rdname performQuery
setMethod("performQuery", "eutil", function(x, method = "GET", ...) {
  method <- match.arg(method, c("GET", "POST"))
  x$perform_query(method = method, ...)
  return(invisible(x))
})

#' database
#' 
#' Retrieve the target database name from an \code{\linkS4class{eutil}} object.
#' 
#' @param x An \code{\linkS4class{eutil}} object.
#' @param ... Further arguments passed on to methods.
#' @return A character string.
#' @seealso
#'    \code{\link{einfo}}, \code{\link{esearch}}, \code{\link{esummary}},
#'    \code{\link{efetch}}, \code{\link{elink}}, \code{\link{epost}},
#'    \code{\link{egquery}}, \code{\link{espell}}, \code{\link{ecitmatch}}.
#' @export
#' @examples
#' \dontrun{
#' e <- esearch("Mus musculus", "taxonomy")
#' database(e)
#' }
setGeneric("database", function(x, ...) standardGeneric("database"))
#' @describeIn database
setMethod("database", "eutil", function(x, ...) x$database())

#' retmode
#' 
#' Get the \dQuote{retrieval mode} of an \code{\linkS4class{eutil}} object
#' It is usually one of \code{xml}, \code{json}, \code{text}, or \code{asn.1}. 
#' It is set to \code{NULL} if \dQuote{retrieval mode} is not supported by an
#' E-Utility.
#' 
#' @param x An \code{\linkS4class{eutil}} object.
#' @param ... Further arguments passed on to methods.
#' @return A character string or \code{NULL}. 
#' @seealso
#'    \code{\link{einfo}}, \code{\link{esearch}}, \code{\link{esummary}},
#'    \code{\link{efetch}}, \code{\link{elink}}, \code{\link{epost}},
#'    \code{\link{egquery}}, \code{\link{espell}}, \code{\link{ecitmatch}}.
#' @export
#' @examples
#' \dontrun{
#' e <- efetch("10090", "taxonomy")
#' retmode(e)
#' }
setGeneric("retmode", function(x, ...) standardGeneric("retmode"))
#' @describeIn retmode
setMethod("retmode", "eutil", function(x, ...) x$retmode())

#' rettype
#' 
#' Get the \dQuote{retrieval type} of an \code{\linkS4class{eutil}} object. See 
#' \href{http://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly}{here}
#' for the available retrieval types for different NCBI databases.
#' 
#' @param x An \code{\linkS4class{eutil}} object.
#' @param ... Further arguments passed on to methods.
#' @return A character string.
#' @seealso
#'    \code{\link{einfo}}, \code{\link{esearch}}, \code{\link{esummary}},
#'    \code{\link{efetch}}, \code{\link{elink}}, \code{\link{epost}},
#'    \code{\link{egquery}}, \code{\link{espell}}, \code{\link{ecitmatch}}.
#' @export
#' @examples
#' \dontrun{
#' e <- esearch("Mus musculus", "taxonomy")
#' rettype(e)
#' }
setGeneric("rettype", function(x, ...) standardGeneric("rettype"))
#' @describeIn rettype
setMethod("rettype", "eutil", function(x, ...) x$rettype())

#' uid
#' 
#' Retrieve the list of UIDs returned by a call to ESearch or ELink.
#' 
#' @param x An \code{\linkS4class{eutil}} object.
#' @param ... Further arguments passed on to methods.
#' @return A character vector.
#' @seealso
#'    \code{\link{esearch}}, \code{\link{elink}}.
#' @export
#' @examples
#' \dontrun{
#' e <- esearch("Mus musculus", "taxonomy")
#' uid(e)
#' }
setGeneric("uid", function(x, ...) standardGeneric("uid"))

#' webenv
#' 
#' Retrieve the Web environment string returned from an ESearch, EPost or ELink call.
#' \code{NA} if the History server was not used.
#' 
#' @param x An \code{\linkS4class{eutil}} object.
#' @param ... Further arguments passed on to methods.
#' @return A character string or \code{NA}.
#' @seealso
#'    \code{\link{einfo}}, \code{\link{esearch}}, \code{\link{esummary}},
#'    \code{\link{efetch}}, \code{\link{elink}}, \code{\link{epost}},
#'    \code{\link{egquery}}, \code{\link{espell}}, \code{\link{ecitmatch}}.
#' @export
#' @examples
#' \dontrun{
#' e <- esearch("Mus musculus", "taxonomy", usehistory = TRUE)
#' webenv(e)
#' }
setGeneric("webenv", function(x, ...) standardGeneric("webenv"))

#' querykey
#' 
#' An integer query key returned by an ESearch, EPost or ELink call if
#' the History server was used. Otherwise \code{NA}.
#' 
#' @param x An \code{\linkS4class{eutil}} object.
#' @param ... Further arguments passed on to methods.
#' @return An integer or \code{NA}.
#' @seealso
#'    \code{\link{einfo}}, \code{\link{esearch}}, \code{\link{esummary}},
#'    \code{\link{efetch}}, \code{\link{elink}}, \code{\link{epost}},
#'    \code{\link{egquery}}, \code{\link{espell}}, \code{\link{ecitmatch}}.
#' @export
#' @examples
#' \dontrun{
#' e <- esearch("Mus musculus", "taxonomy", usehistory = TRUE)
#' querykey(e)
#' }
setGeneric("querykey", function(x, ...) standardGeneric("querykey"))

