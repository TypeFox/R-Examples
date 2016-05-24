#' Search by Vertnet occurrence ID
#'
#' @export
#' @param ids (character) VertNet IDs, one or more. Required.
#' @param compact (logical) Return a compact data frame. That is, remove
#' empty columns. Default: \code{TRUE}
#' @param verbose (logical) Print progress and information messages. Default: \code{TRUE}
#' @param ... Curl arguments passed on to \code{\link[httr]{GET}}
#'
#' @return A list, with data frame of search results, and list of metadata
#' @references \url{http://bit.ly/vertnet-wiki}
#' @details VertNet IDs can be a variety of things, some URIs (i.e., with http://...), 
#' while others start with \code{urn}. 
#'
#' @examples \dontrun{
#' vert_id(ids = "urn:catalog:AUM:Fish:13271")
#' ids <- c("http://arctos.database.museum/guid/MSB:Mamm:56979?seid=1643089", 
#'          "urn:catalog:CM:Herps:116520",
#'          "urn:catalog:AUM:Fish:13271")
#' res <- vert_id(ids)
#' res$data$occurrenceid
#' identical(sort(res$data$occurrenceid), sort(ids))
#' 
#' # Lots
#' out <- vertsearch(taxon = "aves", state = "california", limit = 20)
#' ids <- out$data$occurrenceid
#' res <- vert_id(ids)
#' identical(sort(res$data$occurrenceid), sort(ids))
#' }

vert_id <- function(ids, compact = TRUE, verbose = TRUE, ...) {
  tt <- GET(vurl(), query = list(q = make_id_q(ids)), ...)
  stop_for_status(tt)
  txt <- content(tt, "text", encoding = "UTF-8")
  out <- jsonlite::fromJSON(txt)
  avail <- out$matching_records
  result <- out$recs
  df <- if (NROW(result) == 0) data.frame(NULL, stringsAsFactors = FALSE) else result
  names(df) <- tolower(names(df))
  res <- get_terms()
  df <- merge(res$fullr, df, all = TRUE)[, tolower(res$termlist[,1]) ]
  df <- df[ -NROW(df), ]
  if (compact) { 
    df <- df[ , !sapply(df, function(x) all(is.na(x))) ]
  }
  mssg(verbose, paste("\nLast Query URL: \"", tt$url, "\"", sep = ""))
  mssg(verbose, paste("\nMatching records:", NROW(df), "returned,", avail, "available", sep = " "))
  list(meta = make_meta(out), data = tbl_df(df))
}

make_id_q <- function(x) {
  x <- paste0(sprintf('\"%s\"', unname(unlist(x))), collapse = " OR ")
  x <- list(q = sprintf("iptrecordid:%s", x))
  jsonlite::toJSON(x, auto_unbox = TRUE)
}
