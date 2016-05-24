#' @include eutil.R
NULL

#' @export
.egquery <- setRefClass(
  Class    = "egquery",
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
        methods::show(get_content("xml"))
        tail <- sprintf("EGQuery query for the search term %s.",
                        sQuote(xmlValue("/Result/Term")))
        cat(tail, sep = "\n")
      } else {
        methods::show(get_error())
      }
    }
  )
)

#' \code{egquery} retrieves the number of records in all Entrez databases for
#' a single text query.
#' 
#' @details
#' See the official online documentation for NCBI's
#' \href{http://www.ncbi.nlm.nih.gov/books/NBK25499//#chapter4.EGQuery}{EUtilities}
#' for additional information.
#' 
#' @title egquery - performing a global Entrez search
#' @param term A valid Entrez text query.
#' @return An \code{\linkS4class{egquery}} object.
#' @export
#' @examples
#' ## Determine the number of records for mouse in Entrez.
#' e <- egquery("mouse[orgn]")
#' e
egquery <- function(term) {
  if (missing(term)) {
    stop("No search term provided")
  }
  if (length(term) > 1L) {
    term <- paste(term, collapse=" OR ")
  }
  .egquery('GET', term = term, retmode = 'xml')
}

#' @describeIn content
setMethod("content", "egquery", function(x, as = NULL) {
  callNextMethod(x = x, as = as)
})
