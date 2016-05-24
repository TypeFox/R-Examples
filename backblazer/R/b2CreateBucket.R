#' Create B2 Bucket.
#'
#' \code{b2CreateBucket} creates a new bucket in the user's account on the
#' Backblaze B2 cloud storage product.
#'
#' This function creates a new bucket within the user's account on the Backblaze
#' B2 cloud storage product. Backblaze B2 does not support tree based folder
#' structures as such, meaning all uploaded data is stored in a flat file
#' structure. However, buckets may be created at the top level to aid in content
#' organisation. Further details regarding this API call are available here:
#'
#' \url{https://www.backblaze.com/b2/docs/b2_create_bucket.html}
#'
#' \code{bucketName} and \code{bucketType} are mandatory and must be user
#' defined.
#'
#' @param bucketName Bucket names must be globally unique. No two users may have
#'   buckets named the same. Bucket names may not start with \emph{b2}. Bucket
#'   names must be a minimum of 6 and maximum of 50 characters long. Bucket
#'   names may consist only of letters, numbers and hyphens. Special characters
#'   are invalid.
#' @param bucketType Supported bucket types are \emph{allPublic} and
#'   \emph{allPrivate}.
#' @return If successful a list containing the \code{accountId},
#'   \code{bucketName} and \code{bucketType} will all be echoed back to the
#'   user. Also included in the return list will be a unique \code{bucketId}.
#'
#' @examples
#' \dontrun{
#' b2CreateBucket(bucketName = "this-is-a-uniquely-named-bucket",
#' bucketType = "allPublic")
#' }
#'
#' @export

b2CreateBucket <- function(bucketName, bucketType) {
  # Read Account Authorisation file
  accountAuthorization <- NULL
  accountAuthorization <- readRDS("accountAuthorization.rds")

  # Function options from input, make a dataframe
  accountId <- as.character(accountAuthorization$accountId)
  accountId <- as.data.frame(accountId, stringsAsFactors = FALSE)
  bucketName <- as.data.frame(bucketName, stringsAsFactors = FALSE)
  bucketType <- as.data.frame(bucketType, stringsAsFactors = FALSE)

  # Bind function option data frames together
  vars <- cbind(accountId, bucketName, bucketType)

  # API call
  b2Return <-
    httr::POST(
      paste(
        accountAuthorization$apiUrl,"/b2api/v1/b2_create_bucket", sep = ""
      ), body = jsonlite::toJSON(jsonlite::unbox(vars), pretty = TRUE), httr::add_headers(
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
    jsonlite::fromJSON(httr::content(b2Return, type = "text"))
  }
}
