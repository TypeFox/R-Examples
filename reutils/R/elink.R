#' @include eutil.R
NULL

#' @export
.elink <- setRefClass(
  Class    = "elink",
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

parse_linkset <- function(object) {
  x <- object$get_content("xml")
  DbFrom <- xvalue(x, "/eLinkResult/LinkSet/DbFrom")
  IdList <- xvalue(x, "/eLinkResult/LinkSet/IdList/Id")
  LinkSetDb <- xset(x, "/eLinkResult/LinkSet/LinkSetDb")
  if (length(LinkSetDb) < 1L) {
    lset <- list(
      structure(NA_character_, score = NA_integer_,
                database = object$params[["db"]],
                linkName = paste0(DbFrom, "_", object$params[["db"]]),
                class = c("entrez_link", "character"))
      )
  } else {
    lset <- lapply(LinkSetDb, function(lsd) {
      lsd <- XML::xmlDoc(lsd)
      uid <- xvalue(lsd, "/LinkSetDb/Link/Id")
      score <- xvalue(lsd, "/LinkSetDb/Link/Score", as = "integer")
      dbTo <- xvalue(lsd, "/LinkSetDb/DbTo") %|char|% NA_character_
      linkName <- xvalue(lsd, "/LinkSetDb/LinkName") %|char|% NA_character_
      XML::free(lsd)
      structure(uid, score = score, database = dbTo, linkName = linkName,
                class = c("entrez_link", "character"))
    })
    lset <- compactNA(lset)
  }

  lnm <- vapply(lset, attr, "linkName", FUN.VALUE = "")
  x <- structure(
    lset,
    names = lnm,
    database = DbFrom,
    uid = IdList,
    class = c("entrez_linkset", "list")
  )
  x
}

#' Class \code{"entrez_linkset"}
#'
#' A list containing a set of links as returned by a call to \code{\link{elink}}.
#' Each element of the list is a character vector of UIDs of class 
#' \code{"entrez_link"} with three attributes:
#' \describe{
#'    \item{\code{score}:}{Similarity scores between query UIDs and the linked UIDs}
#'    \item{\code{database}:}{The destination database of the ELink query.}
#'    \item{\code{linkName}:}{Name of the retrieved Entrez link of the form
#'    \emph{dbFrom_dbTo_subset}}
#' }
#' 
#' An \code{"entrez_linkset"} has two global attributes:
#' \describe{
#'    \item{\code{uid}:}{The input UIDs.}
#'    \item{\code{database}:}{The database containing the input UIDs.}
#' }
#' 
#' @keywords classes internal
#' @name entrez_linkset-class
#' @examples
#' ###
setOldClass("entrez_linkset")

#' @export
print.entrez_link <- function(x, ...) {
  db <- strsplit(attr(x, "linkName"), "_")[[1]]
  dbFrom <- db[1]
  dbTo <- paste0(db[-1], collapse = "_")
  cat(sprintf("List of linked UIDs from database %s to %s.\n", sQuote(dbFrom), sQuote(dbTo)))
  print(format(x))
  invisible()
}

#' @export
"[.entrez_link" <- function(x, i, j, ..., drop=TRUE) {
  out <- NextMethod(...)
  attr(out, "score") <- attr(x, "score")[i]
  attr(out, "database") <- attr(x, "database")
  attr(out, "linkName") <- attr(x, "linkName")
  class(out) <- c("entrez_link", "character")
  out
}

#' @describeIn database
setMethod("database", "entrez_linkset", function(x, ...) attr(x, "database"))

#' @describeIn uid
setMethod("uid", "entrez_linkset", function(x, ...) attr(x, "uid"))

#' linkset
#' 
#' Retrieve a linkset from an \code{\linkS4class{elink}} object.
#' 
#' @param x An \code{\linkS4class{elink}} object.
#' @param linkname (optional) Name of the Entrez link to retrieve. Every link in
#' Entrez is given a name of the form \emph{dbFrom_dbTo_subset}. If \code{NULL},
#' all available links are retrieved from the object.
#' @param ... Further arguments passed on to methods.
#' @return A list.
#' @export
#' @examples
#' \dontrun{
#' ## Find related articles to PMID 20210808 and xtract linked UIDs from the
#' ## "pubmed" to "pubmed_reviews" link
#' x <- elink("20210808", dbFrom = "pubmed", dbTo = "pubmed", cmd = "neighbor_score")
#' linkset(x, "pubmed_pubmed_reviews")
#' }
setGeneric("linkset", function(x, linkname = NULL, ...) standardGeneric("linkset"))
#' @describeIn linkset
setMethod("linkset", "entrez_linkset", function(x, linkname = NULL, ...) {
  if (!is.null(linkname)) {
    ans <- x[linkname]
    if (is.null(ans)) {
      return(NULL)
    }
    if (length(ans) == 1) {
      ans  <- ans[[1]]
    }
    ans
  } else {
    attr(x, "database") <- NULL
    attr(x, "uid") <- NULL
    unclass(x)
  }
})

#' @export
print.entrez_linkset <- function(x, ...) {
  db <- database(x)
  dbTo <- vapply(x, attr, 'database', FUN.VALUE = "")
  cat(sprintf("ELink query from database %s to destination database %s.\n",
              sQuote(db), sQuote(unique(dbTo))))
  cat("Query UIDs:\n")
  print(format(uid(x)))
  cat("Summary of LinkSet:\n")
  lnames <- vapply(x, attr, 'linkName', FUN.VALUE = "")
  llen <- if (all(is.na(x))) {
    0
  } else vapply(x, length, 0)
  print(format(data.frame(DbTo = unname(dbTo), LinkName = unname(lnames),
                          LinkCount = unname(llen))))
}


#' \code{elink} generates a list of UIDs in a specified Entrez database that
#' are linked to a set of input UIDs in either the same or another
#' database. For instance, the ELink utility can find Entrez gene records
#' linked to records in Entrez Protein.
#' 
#' @title elink - finding related data through Entrez links
#' @details
#' See the official online documentation for NCBI's
#' \href{http://www.ncbi.nlm.nih.gov/books/NBK25499//#chapter4.ELink}{EUtilities}
#' for additional information.
#' 
#' If \code{dbTo} and \code{dbFrom} are set to the same database, ELink will
#' return neighbors within that database.
#' 
#' Elink commands (cmd) specify the function that elink will perform.
#' Available commands are:
#' \itemize{
#'   \item{"\strong{neighbor}" }{(Default) ELink returns a set of UIDs in dbTo
#'   linked to the input UIDs in dbFrom.}
#'   \item{"\strong{neighbor_score}" }{ELink returns a set of UIDs within the
#'   same database as the input UIDs along with similarity scores.}
#'   \item{"\strong{neighbor_history}" }{ELink posts the output UIDs to the
#'   Entrez History server and returns a query_key and WebEnv parameter.
#'   Alternatively this is achieved by setting \code{usehistory=TRUE}}
# '   \item{"\strong{acheck}" }{ELink lists all links available for a set of UIDs.}
# '   \item{"\strong{ncheck}" }{ELink checks for the existence of links
# '   \emph{within the same database} for a set of UIDs.}
# '   \item{"\strong{lcheck}" }{Elink checks for the existence of external links
# '   (LinkOuts) for a set of UIDs.}
# '   \item{"\strong{llinks}" }{For each input UID, ELink lists the URLs and
# '   attributes for the LinkOut providers that are not libraries.}
# '   \item{"\strong{llinkslib}" }{For each input UID, ELink lists the URLs and
# '   attributes for all LinkOut providers including libraries.}
# '   \item{"\strong{prlinks}" }{ELink lists the primary LinkOut provider for
# '   each input UID.}
#' }
#' 
#' @param uid (Required) A character vector of UIDs.
#' @param dbFrom Initial database containing the UIDs in the input list.
#' @param dbTo Destination database from which to retrieve linked UIDs. If
#' not provided links will be sought in the database containing the input UIDs.
#' @param linkname Name of the Entrez link to retrieve. Every link in
#' Entrez is given a name of the form \emph{dbFrom_dbTo_subset}.
#' @param usehistory If \code{TRUE} search results are stored directly in
#' the user's Web environment so that they can be used in subsequents 
#' calls to \code{\link{esummary}} or \code{\link{efetch}}.
#' @param cmd ELink command mode (default: 'neighbor'). See Details.
#' @param correspondence if \code{TRUE} correspondence between query UIDs and
#' destination UIDs is preserved.
#' @param querykey Query key.
#' @param webenv Web Environment.
#' @param term Search query to limit the output set of linked UIDs.
#' @param holding Name of LinkOut provider.
#' @param datetype Type of date to limit the search. One of 'mdat'
#' (modification date), 'pdat' (publication date) or 'edat' (Entrez date).
#' @param reldate umber of days back for which search items are
#' returned.
#' @param mindate Minimum date of search range. Format YYYY/MM/DD.
#' @param maxdate Maximum date of search range. Format YYYY/MM/DD.
#' @return An \code{\linkS4class{elink}} object.
#' @export
#' @seealso
#' Combine calls to ELink with other EUtils:
#' \code{\link{esummary}}, \code{\link{efetch}}.
#' @seealso
#' Accessor methods:
#' \code{\link{content}}, \code{\link{getUrl}}, \code{\link{getError}},
#' \code{\link{database}}, \code{\link{uid}}, \code{\link{linkset}}, 
#' @examples
#' ## Find one set of Gene IDs linked to nuccore GIs 34577062 and 24475906
#' e <- elink(c("927442695", "312836839"), dbFrom = "nuccore", dbTo = "gene")
#' e
#' 
#' \dontrun{
#' ## Find related articles to PMID 20210808
#' p <- elink("20210808", dbFrom = "pubmed", dbTo = "pubmed")
#' p
#' 
#' ## Extract linked UIDs from the "pubmed" to "pubmed_reviews" link
#' linkset(p, "pubmed_pubmed_reviews")
#' 
#' ## or
#' p["pubmed_pubmed_reviews"]
#' 
#' ## retrive the abstracts for the first five linked reviews
#' abstracts <- efetch(p["pubmed_pubmed_reviews"][1:5], rettype = "abstract")
#' }
elink <- function(uid, dbFrom = NULL, dbTo = NULL, linkname = NULL,
                  usehistory = FALSE, cmd = "neighbor",
                  correspondence = FALSE, querykey = NULL, webenv = NULL,
                  term = NULL, holding = NULL, datetype = NULL,
                  reldate = NULL, mindate = NULL, maxdate = NULL) {
  ## extract query parameters
  params <- parse_params(uid, dbFrom, querykey, webenv)
  
  ## set dbTo=dbFrom if no dbTo is provided
  if (is.null(dbTo) && !grepl(pattern = "check$|links", cmd)) {
    dbTo <- params$db
  }
  if (usehistory) {
    cmd <- "neighbor_history"
  }
  if  (correspondence && !is.null(params$uid)) {
    uid  <- paste0(params$uid, collapse = "&id=")
  } else {
    uid <- .collapse(params$uid)
  }
  .elink(method = if (length(params$uid) < 100) "GET" else "POST",
         id = uid, db = dbTo, dbFrom = params$db, cmd = cmd, query_key = params$querykey,
         WebEnv = params$webenv, linkname = linkname, term = term, holding = holding,
         datetype = datetype, reldate = reldate, mindate = mindate, maxdate = maxdate,
         retmode = 'xml')
}

#' @describeIn content
setMethod("content", "elink", function(x, as = NULL) {
  callNextMethod(x = x, as = as)
})

#' ELink Accessors
#' 
#' Extract UIDs from an \code{\link{elink}} object.
#'
#' @param x An \code{\linkS4class{elink}} object.
#' @param i Numeric or character indices.
#' @param j Ignored.
#' @return A \code{\linkS4class{entrez_linkset}} object.
#' @rdname extract-elink
#' @examples
#' \dontrun{
#' e <- elink(c("34577062", "24475906"), dbFrom = "nuccore")
#' e[1]
#' }
setMethod("[", c(x = "elink", i = "ANY", j = "missing"), function(x, i, j) {
  linkset(x, i)
})

#' @rdname extract-elink
setMethod("[", c("elink", "character"), function(x, i) {
  linkset(x, i)
})

#' @describeIn linkset
setMethod("linkset", "elink", function(x, linkname = NULL, ...) {
  linkset(x$get_content("parsed"), linkname = linkname, ...)
})

#' @describeIn uid
setMethod("uid", "elink", function(x, ...) {
  uid(x$get_content("parsed"))
})

