#' @include eutil.R
NULL

#' @export
.einfo <- setRefClass(
  Class    = "einfo",
  contains = "eutil",
  methods  = list(
    initialize = function(method = "GET", ...) {
      callSuper()
      perform_query(method = method, ...)
      if (no_errors()) {
        errors$check_errors(.self)
      }
    },
    show_xml = function() {
      if (is.null(database())) {
        cat("List of Entrez databases\n")
        methods::show(xmlValue("/eInfoResult/DbList/DbName"))
      } else {
        cat(sprintf("Overview over the Entrez database %s.\n",
                    sQuote(xmlValue("/eInfoResult/DbInfo/MenuName"))))
        x <- get_content("parsed")
        nm <- names(x)
        fnm <- ellipsize(paste0(names(x$Fields), collapse = "; "), offset = 10)
        lnm <- ellipsize(paste0(names(x$Links), collapse = "; "), offset = 9)
        showme <- sprintf("  %s: %s", nm, c(vapply(x[1:6], as.character, ""),
                                            fnm, lnm))
        cat(showme, sep = "\n")
      }
    },
    show_json = function() {
      methods::show(get_content("json"))
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

parse_einfo <- function(object) {
  retmode <- object$retmode()
  if (retmode == "xml") {
    x <- object$get_content("xml")
    if (is.null(object$database())) {
      xvalue(x, "/eInfoResult/DbList/DbName")
    } else {
      list(
        dbName      = xvalue(x, "/eInfoResult/DbInfo/DbName"),
        MenuName    = xvalue(x, "/eInfoResult/DbInfo/MenuName"),
        Description = xvalue(x, "/eInfoResult/DbInfo/Description"),
        DbBuild     = xvalue(x, "/eInfoResult/DbInfo/DbBuild"),
        Count       = xvalue(x, "/eInfoResult/DbInfo/Count", "integer"),
        LastUpdate  = xvalue(x, "/eInfoResult/DbInfo/LastUpdate", "POSIXlt"),
        Fields      = extract_df(x, "/eInfoResult/DbInfo/FieldList/Field/"),
        Links       = extract_df(x, "/eInfoResult/DbInfo/LinkList/Link/")
      )
    }
  } else if (retmode == "json") {
    jsonlite::fromJSON(object$get_content("json"))
  }
}

extract_df <- function(x, path) {
  nm <- unique(xname(x, paste0(path, "child::node()")))        
  if (!all(is.na(nm))) {
    nodes <- xset(x, paste0(path, "*"))
    list <- split(vapply(nodes, XML::xmlValue, ""), nm)
    data.frame(stringsAsFactors = FALSE, list)[, nm]
  } else {
    data.frame()
  }
}

#' \code{einfo} queries the NCBI EInfo utility to retrieve the names of all
#' valid Entrez databases, or, if \code{db} is provided, to retrieve statistics
#' for a single database, including lists of indexing fields and available link
#' names. Version 2.0 data is requested by default. 
#' 
#' @details
#' See the official online documentation for NCBI's
#' \href{http://www.ncbi.nlm.nih.gov/books/NBK25499//#chapter4.EInfo}{EUtilities}
#' for additional information.
#' 
#' @title einfo - getting database statistics and search fields
#' @param db A valid NCBI database name. If \code{NULL}, a list of all current
#' NCBI databases is returned.
#' @param version Specifies version 2.0 EInfo XML. Set to \code{NULL} for the
#' older version.
#' @param retmode 'xml' (default) or 'json'.
#' @return An \code{\linkS4class{einfo}} object.
#' @seealso
#' \code{\link{content}}, \code{\link{getUrl}}, \code{\link{getError}}.
#' @export
#' @examples
#' \dontrun{
#' ## Fetch a list of all current Entrez database names
#' einfo()
#' 
#' ## Fetch statistics for an Entrez database and parse
#' ## the data into a data.frame
#' x <- einfo("gene")
#' if (x$no_errors()) {
#'   content(x, "parsed")
#' }
#' 
#' 
#' 
#' ## Fetch statistics for an Entrez database in JSON format
#' ## and parse the data into a list
#' x <- einfo("pubmed", retmode = "json")
#' if (x$no_errors()) {
#'   content(x, "parsed")
#' }
#' }
einfo <- function(db = NULL, version = "2.0", retmode = "xml") {
  retmode <- match.arg(retmode, c("xml", "json"))
  assertthat::assert_that(is.null(db) || assertthat::is.string(db))
  if (!is.null(db)) {
    assertthat::assert_that(is.null(version) || version == "2.0")
  } else {
    version <- NULL
  }
  
  .einfo(method = "GET", db = db, version = version, retmode = retmode)
}

#' @describeIn content
setMethod("content", "einfo", function(x, as = NULL) {
  callNextMethod(x = x, as = as)
})

#' EInfo accessors
#' 
#' Extract parts of a parsed \code{\linkS4class{einfo}} object.
#' 
#' @param x An \code{\linkS4class{einfo}} object.
#' @param i Numeric or character indices specifying the elements to extract.
#' @param j Ignored.
#' @return A list.
#' @seealso \code{\link[base]{Extract}}
#' @rdname extract-einfo
#' @examples
#' \dontrun{
#' e <- einfo("pubmed")
#' e[1:5]
#' e["Description"]
#' e[["Links"]]
#' 
#' e2 <- einfo("pubmed", retmode = 'json')
#' e2[["header"]]
#' e2[["einforesult"]][["dbinfo"]][["description"]]
#' }
setMethod("[", c(x = "einfo", i = "ANY", j = "missing"), function(x, i, j) {
  content(x, as = "parsed")[i]
})

#' @rdname extract-einfo
setMethod("[[", "einfo", function(x, i) {
  content(x, as = "parsed")[[i]]
})
