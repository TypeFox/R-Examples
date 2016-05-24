#' Search Rodham's emails
#'
#' @description Search Hillary Rodham Clinton's \emph{personal} emails.
#'
#' @param subject Filter by subject, defaults to \code{NULL}(no filter). If
#' \code{internal = TRUE} then matches pattern, if \code{internal = FALSE} then
#' looks for exact match.
#' @param to Filter by Receiver, defaults to \code{NULL}(no filter).
#' @param from Filter by Sender, defaults to \code{NULL}(no filter).
#' @param start Filter by date range, defaults to \code{NULL}(no filter).
#' @param end Filter by date range, defaults to \code{NULL}(no filter).
#' @param internal if \code{TRUE} (default) searches the internal data set
#' (see \code{data(emails)}), if \code{FALSE} fetches the data through
#' the Wall Street journal API. \code{data(emails)} is equivalent to internal
#' \code{TRUE}
#'
#' @details There are a total of 29444 emails ranging from \code{2009-08-14} to
#' \code{2014-08-13}, please consider leaving internal to \code{TRUE} to not
#' hammer the Wall Street Journal's API. \code{internal = TRUE} is equivalent
#' to \code{\link{emails}}.
#'
#' @examples
#' \dontrun{
#' emails <- search_emails()
#'
#' # only emails on cuba
#' emails <- search_emails(subject = "Cuba")
#'
#' # only emails from Jake Sullivan since 2014
#' j_s <- search_emails(from = "Jake Sullivan", start = as.Date("2014-01-01"))
#' }
#'
#' @author John Coene \email{jcoenep@@gmail.com}
#'
#' @export
search_emails <- function(subject = NULL, to = NULL, from = NULL, start = NULL,
                         end = NULL, internal = TRUE){
  if (internal == FALSE) { # encode subject if fetching from WSJ
    if (!is.null(subject)) {
      subject <- URLencode(toupper(subject))
    }
    uri <- paste0("http://graphics.wsj.com/hillary-clinton-email-documents/api/",
                  "search.php?subject=", subject,
                  "&text=&to=", to,
                  "&from=", from,
                  "&start=", start,
                  "&end=", end,
                  "&sort=docDate&order=desc&docid=&limit=30000&offset=0")
    json <- jsonlite::fromJSON(uri)
    if (json$total > 0) {
      emails <- json$rows
      cat(json$total, "emails returned")
    } else if (json$total == 0){
      cat("No emails returned")
      emails <- data.frame()
    }
  } else {
    if (!is.null(subject)) {
      emails <- emails[grep(toupper(subject), emails$subject),]
    }
    if (!is.null(to)) {
      emails <- emails[emails$to == to,]
    }
    if (!is.null(from)) {
      emails <- emails[emails$from == from,]
    }
    if (!is.null(start)) {
      emails <- emails[emails$docDate <= as.Date(start),]
    }
    if (!is.null(end)) {
      emails <- emails[emails$docDate >= as.Date(end),]
    }
  }
  return(emails)
}
