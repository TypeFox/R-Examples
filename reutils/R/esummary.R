#' @include eutil.R
#' @include parse-params.R
NULL

#' @export
.esummary <- setRefClass(
  Class    = "esummary",
  contains = "eutil",
  methods  = list(
    initialize = function(method, ...) {
      callSuper()
      perform_query(method = method, ...)
      if (no_errors()) {
        errors$check_errors(.self)
      }
    },
    show_xml = function() {
      methods::show(get_content("xml"))
      tail <- sprintf("ESummary query using the database %s.", sQuote(database()))
      cat(tail, sep="\n")
    },
    show_json = function() {
      methods::show(get_content("json"))
      tail <- sprintf("ESummary query using the database %s.", sQuote(database()))
      cat(tail, sep="\n")
    },
    show = function() {
      cat("Object of class", sQuote(eutil()), "\n")
      if (no_errors()) {
        switch(retmode(), xml = show_xml(), json = show_json())
      } else {
        methods::show(get_error())
      }
    } 
  )
)

#' \code{esummary} performs calls to the NCBI ESummary utility to retrieve document
#' summaries (DocSums) for a list of primary UIDs or for a set of UIDs stored in the
#' user's web environment (using the Entrez History server).
#' 
#' @details
#' See the official online documentation for NCBI's
#' \href{http://www.ncbi.nlm.nih.gov/books/NBK25499//#chapter4.ESummary}{EUtilities}
#' for additional information.
#'
#' @title esummary - downloading Document Summaries
#' @param uid (Required)
#' List of UIDs provided either as a character vector, as an
#' \code{esearch} or \code{elink} object, or by reference to a Web
#' Environment and a query key obtained directly from objects returned
#' by previous calls to \code{\link{esearch}}, \code{\link{epost}} or
#' \code{\link{elink}}.
#' If UIDs are provided as a plain character vector, \code{db} must be
#' specified explicitly, and all of the UIDs must be from the database
#' specified by \code{db}.
#' @param db (Required only when \code{id} is a character vector of UIDs)
#' Database from which to retrieve DocSums.
#' @param retstart Numeric index of the first DocSum to be retrieved
#' (default: 1).
#' @param retmax Total number of DocSums from the input set to be retrieved
#' (maximum: 10,000).
#' @param querykey An integer specifying which of the UID lists attached
#' to a user's Web Environment will be used as input to \code{efetch}.
#' (Usually obtained drectely from objects returned by previous
#' \code{\link{esearch}}, \code{\link{epost}} or \code{\link{elink}} calls.)
#' @param webenv A character string specifying the Web Environment that
#' contains the UID list. (Usually obtained directely from objects returned
#' by previous \code{\link{esearch}}, \code{\link{epost}} or
#' \code{\link{elink}} calls.)
#' @param retmode Retrieval mode. (default: 'xml', alternative: 'json')
#' @param version If "2.0" \code{esummary} will retrieve version 2.0
#' ESummary XML output.
#' @return An \code{\linkS4class{esummary}} object.
#' @seealso
#' \code{\link{content}}, \code{\link{getUrl}}, \code{\link{getError}},
#' \code{\link{database}}.
#' @export
#' @examples
#' ## Retrieve the Document Summary information for a set of
#' ## UIDs frome the Gene datanase.
#' ds <- esummary(c("828392", "790", "470338"), "gene")
#' ds
#' 
#' \dontrun{
#' ## parse the XML into a data frame
#' df <- content(ds, "parsed")
#' df
#' 
#' ## use XPath expressions to extract nodes of interest
#' ds['//TaxID']
#' }
esummary <- function(uid, db = NULL, retstart = 1, retmax = 10000,
                      querykey = NULL, webenv = NULL, retmode = 'xml',
                     version = "2.0") {
  ## extract query parameters
  params <- parse_params(uid, db, querykey, webenv)
  retmode <- match.arg(retmode, c('xml', 'json'))
  if (retmax > 10000) {
    stop("Number of DocSums to be downloaded should not exceed 10,000.", call.=FALSE)
  }
  .esummary(method = if (length(params$uid) < 100) "GET" else "POST",
            db = params$db, id = .collapse(params$uid),
            query_key = params$querykey, WebEnv = params$webenv, 
            retstart = retstart, retmax = retmax, retmode = retmode,
            version = if (version == "2.0") "2.0" else NULL)
}

#' @describeIn content
setMethod("content", "esummary", function(x, as = NULL) {
  callNextMethod(x = x, as = as)
})

#' ESummary accessors
#' 
#' Extract XML nodes from an \code{\linkS4class{esummary}} object.
#' 
#' @param x An \code{\linkS4class{esummary}} object.
#' @param i An XPath expression.
#' @return An XML node set.
#' @rdname sub-esummary
#' @export
#' @examples
#' \dontrun{
#' ds <- esummary("470338", "protein")
#' ds["//Slen/node()"]
#' 
#' library("XML")
#' as.numeric(xmlValue(ds[["//Slen"]]))
#' }
setMethod("[", c("esummary", "character"), function(x, i) {
  x$xmlSet(i)  
})

#' @rdname sub-esummary
#' @export
setMethod("[[", c("esummary", "character"), function(x, i) {
  ans <- x[i]
  if (length(ans) > 1) {
    warning(length(ans), " elements in node set. Returning just the first!")
  }
  ans[[1]]
})
