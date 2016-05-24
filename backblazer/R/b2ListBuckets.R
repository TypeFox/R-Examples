#' List B2 Buckets.
#'
#' \code{b2ListBuckets} lists buckets associated with an account, in
#' alphabetical order by bucket ID.
#'
#' This function lists buckets associated with a user's account on the Backblaze
#' B2 cloud storage product, in alphabetical order by bucket ID. Further details
#' regarding this API call are available here:
#'
#' \url{https://www.backblaze.com/b2/docs/b2_list_buckets.html}
#'
#' \code{b2ListBuckets} does not support any function options. However, the
#' function must first be authorised through the \code{b2AuthorizeAccount}
#' function.
#'
#' @return If successful a data frame will be returned containing
#'   \code{accountId}, \code{bucketId}, \code{bucketName} and \code{bucketType}
#'   for all existing buckets in the specific user's account.
#'
#' @examples
#' \dontrun{
#' b2ListBuckets()
#' }
#'
#' @export

b2ListBuckets <- function() {
  # Read Account Authorisation file
  accountAuthorization <- NULL
  accountAuthorization <- readRDS("accountAuthorization.rds")

  # API call
  b2Return <-
    httr::GET(
      paste(
        accountAuthorization$apiUrl,"/b2api/v1/b2_list_buckets?accountId=",accountAuthorization$accountId,sep =
          ""
      ), httr::add_headers(
        'Authorization' = as.character(accountAuthorization$authorizationToken)
      )
    )

  # Check for bad authorisation and sent message
  if (httr::status_code(b2Return) != "200") {
    badReturn <- jsonlite::fromJSON(httr::content(b2Return,type = "text"))
    stop(
      "\nSomething went wrong. Please check the function options to ensure valid values. \n",
      "\nStatus Code: ", badReturn$code, "\nMessage: ", badReturn$message
    )

  } else {
    # Output as dataframe
    b2Return <- jsonlite::fromJSON(httr::content(b2Return, type = "text"))
    b2Return[[1]]
  }
}
