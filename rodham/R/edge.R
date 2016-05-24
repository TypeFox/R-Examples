#' Network the treacherous
#'
#' @description Builds edge list from emails using \code{from} and \code{to},
#' edge source and target respectively.
#'
#' @param emails Data frame of emails as returned by \code{\link{search_emails}},
#' defaults to \code{emails} (see \code{\link{emails}})
#' @param ... any additional column to keep as meta-data
#'
#' @examples
#' \dontrun{
#' emails <- search_emails()
#'
#' edges <- edges_emails(emails)
#' }
#'
#' @seealso \code{\link{search_emails}}
#'
#' @author John Coene \email{jcoenep@@gmail.com}
#'
#' @export
edges_emails <- function(emails = emails, ...){
  if (missing(emails)) {
    stop("Missing emails, see search_emails")
  }
  emails <- emails[emails$to != "",] # filter
  emails <- emails[emails$from != "",]
  emails$to <- trimws(emails[, "to"]) # trim white space
  emails$from <- trimws(emails[, "from"])
  clean <- emails[with(emails, !grepl(";", to) & !grepl(";", from)),] # split
  raw <- emails[with(emails, grepl(";", to) | grepl(";", from)),]
  tail <- raw2clean(raw) # process raw
  edges <- rbind.data.frame(clean, tail) # bind
  src_tgt <- data.frame(from = edges$from, to = edges$to) # build table
  args <- unlist(list(...))
  if(!is.null(args)){ # if meta-data
    edges$to <- NULL
    edges$from <- NULL
    edges <- as.data.frame(edges)
    src_tgt <- cbind.data.frame(src_tgt,
                                edges[, which(names(edges) %in% args)])
    names(src_tgt)[3:ncol(src_tgt)] <- args
    src_tgt <- plyr::count(src_tgt, c("from", "to", args))
  } else {
    src_tgt <- plyr::count(src_tgt)
  }
  src_tgt <- src_tgt[order(-src_tgt[ncol(src_tgt)]),] # arrange by weight
  return(src_tgt)
}
