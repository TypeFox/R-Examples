#' @include eutil.R
#' @include parse-params.R
NULL

#' @export
.epost <- setRefClass(
  Class    = "epost",
  contains = "eutil",
  methods  = list(
    initialize = function(method, ...) {
      callSuper()
      perform_query(method = method, ...)
      if (errors$all_empty()) {
        errors$check_errors(.self)
      }
    },
    show = function() {
      cat("Object of class", sQuote(eutil()), "\n")
      if (no_errors()) {
        methods::show(get_content("parsed"))
      } else {
        methods::show(get_error())
      }
    }
  )
)

parse_epost <- function(object) {
  if (object$no_errors()) {
    x <- object$get_content("xml")
    structure(
      NA_character_,
      ## Attributes
      retmax   = NA_integer_,
      retstart = NA_integer_,
      count    = count_char(",", object$params$id) + 1,
      query_translation = NA_character_,
      querykey = xvalue(x, '/ePostResult/QueryKey', as = 'numeric'),
      webenv   = xvalue(x, '/ePostResult/WebEnv'),
      database = object$database(),
      class = c("entrez_uid", "character")
    )
  } else {
    structure(NA_character_, database=NA_character_, class=c("entrez_uid", "character"))
  }
}

#' \code{epost} uses the Entrez EPost utility to upload primary UIDs to the Entrez History server
#' or append a list of UIDs to an existing set of UIDs attached to a Web Environment.
#'
#' \code{epost} returns an integer label called a query key and an encoded cookie string
#' called a Web environment. \code{epost} objects can then be used instead of a UID
#' list in subsequent calls to \code{\link{esummary}}, \code{\link{efetch}}, or
#' \code{\link{elink}}.
#' 
#' @details
#' See the official online documentation for NCBI's
#' \href{http://www.ncbi.nlm.nih.gov/books/NBK25499//#chapter4.EPost}{EUtilities}
#' for additional information.
#' 
#' @title epost - uploading UIDs to Entrez
#' @param uid (Required) List of UIDs provided as a character or as an \code{esearch}
#' object.
#' @param db (Required if \code{uid} is a character vector) Database containing the
#' UIDs in the input list.
#' @param webenv (Optional) Web Environment. If provided, this parameter
#' specifies the Web Environment that will receive the UIDs sent by
#' \code{epost}. \code{epost} will create a new query key associated with that
#' Web Environment. The \code{webenv} value is usually returned by a previous
#' call to \code{\link{esearch}}, \code{\link{epost}} or \code{\link{elink}}.
#' If no \code{webenv} parameter is provided, the EPost utility will create a
#' new Web Environment and post the UIDs to query key 1.
#' @return An \code{\linkS4class{epost}} object.
#' @export
#' @examples
#' ## post a list of protein GIs to the Entrez History server
#' gi <- c("194680922", "50978626", "28558982", "9507199", "6678417")
#' p <- epost(gi, "protein")
#' p
epost <- function(uid, db = NULL, webenv = NULL) {
  params <- parse_params(uid, db)
 .epost(method = if (length(params$uid) < 100) "GET" else "POST",
               id = .collapse(params$uid), db = params$db, WebEnv = webenv,
               retmode = 'xml')
}

#' @describeIn content
setMethod("content", "epost", function(x, as = NULL) {
  callNextMethod(x = x, as = as)
})

#' @describeIn webenv
setMethod("webenv", "epost", function(x, ...) webenv(x$get_content("parsed")))

#' @describeIn querykey
setMethod("querykey", "epost", function(x, ...) querykey(x$get_content("parsed")))
