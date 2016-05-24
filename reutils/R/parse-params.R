#' @include utils.R
NULL

.parse_params <- function(uid) {
  if (class(uid)[1] %in% c("epost", "esearch", "elink")) {
    uid <- uid$get_content("parsed")
  }
  if (is(uid, "entrez_linkset") || is(uid, "list")) {
    uid <- merge_linkset(uid)
  }
  params <- list()
  params$uid <- as.character(uid) %|na|% NULL
  params$db <- attr(uid, "database")
  params$querykey <- attr(uid, "querykey") %|na|% NULL
  params$webenv <- attr(uid, "webenv") %|na|% NULL
  if (!is.null(params$webenv)) {
    # if we use the History Server set count to the
    # total number of UIDs stored on the history server
    params$count <- attr(uid, "count")
  } else {
    # set count to the number of UIDs included in the XML output
    # or to the actual number of uids passed into .parse_params
    params$count <- attr(uid, "retmax") %||% length(uid)
  }
  params   
}


parse_params <- function(uid, db = NULL, querykey = NULL, webenv = NULL) {
  ## uids may only be missing if WebEnv and query_key are provided
  if (missing(uid) && (is.null(querykey) || is.null(webenv))) {
    stop("No UIDs provided", call. = FALSE)
  }
  ## if WebEnv and query_key are provided, db must also be provided
  if (!is.null(querykey) && !is.null(webenv) && is.null(db)) {
    stop("No database name provided", call. = FALSE)
  }
  if (missing(uid)) {
    ## if WebEnv and query_key is provided by the user set uid=NULL, count=0, 
    ## restricted to 500.
    params <- list(uid = NULL, db = db, querykey = querykey, webenv = webenv, count = 0)
  } else {
    params <- .parse_params(uid)
    ## abort if uid did not contain db 
    params$db <- db %||% params$db
    if (is.null(params$db)) {
      stop("No database name provided", call. = FALSE)
    }
  }
  
  params
}
