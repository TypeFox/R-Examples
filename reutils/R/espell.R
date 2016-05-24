#' @include eutil.R
NULL

#' @export
.espell <- setRefClass(
  Class    = "espell",
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
      methods::show(get_content("xml"))
    }
  )
)

#' For a text query retrieve an XML containing the original query
#' and spelling suggestions.
#' 
#' @title espell - retrieving spelling suggestions
#' @param term An Entrez text query.
#' @param db An Entrez database.
#' @return An \code{\linkS4class{espell}} object.
#' @export
#' @examples
#' e <- espell("Chlamidia")
#' e
espell <- function(term, db = "pubmed") {
  if (missing(term)) {
    stop("No query term provided", call. = FALSE)
  }
  if (!nzchar(db)) {
    stop("No database provided", call. = FALSE)
  }
  if (length(term) > 1L) {
    term <- paste(term, collapse=" OR ")
  }
  .espell(method = if (nchar(term) < 100) "GET" else "POST",
          term = .escape(term), db = db, retmode = 'xml')
}

#' @describeIn content
setMethod("content", "espell", function(x, as = NULL) {
  callNextMethod(x = x, as = as)
})
