#' Hillary Rodham Clinton emails
#'
#' A dataset containing 29444 emails from/to Hillary Rodham Clinton
#' sent/received between 2009-08-14 and 2014-08-13.
#'
#' @format A data frame with 29444 rows and 9 variables:
#' \describe{
#'   \item{docID}{Primary key}
#'   \item{docDate}{Date when document was sent or received}
#'   \item{to}{Who the emails was sent to}
#'   \item{from}{Who the email is received from}
#'   \item{originalTo}{From whom the email originally comes from}
#'   \item{originalFrom}{To whom the email was originally sent to}
#'   \item{subject}{Subject of the email}
#'   \item{interesting}{Rating, relevancy of email}
#'   \item{not_interesting}{Rating, irrelevancy of email}
#' }
#' @source \url{http://graphics.wsj.com/hillary-clinton-email-documents/api/search.php?subject=&to=&from=&start=&end=&sort=docDate&order=desc&docid=&limit=27159&offset=0}
"emails"
