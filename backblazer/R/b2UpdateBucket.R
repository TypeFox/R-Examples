#' Update B2 Bucket.
#'
#' \code{b2UpdateBucket} modifies the bucket type of an existing bucket in the
#' user's account.
#'
#' This function modifies the bucket type of an existing bucket, within the
#' user's account on the Backblaze B2 cloud storage product. This function can
#' be used to allow everyone to download the contents of the bucket without
#' providing any authorisation, or to prevent anyone from downloading the
#' contents of the bucket without providing a bucket \code{authToken}. Further
#' details regarding this API call are available here:
#'
#' \url{https://www.backblaze.com/b2/docs/b2_update_bucket.html}
#'
#' \code{bucketId} and \code{bucketType} are mandatory and must be user defined.
#'
#' @param bucketId The unique identifier of the bucket to be updated. A list of
#'   all the user's bucket IDs may be found using the \code{b2_list_buckets}
#'   function in this package.
#' @param bucketType Supported bucket types are \emph{allPublic} and
#'   \emph{allPrivate}.
#' @return If successful a list containing the \code{accountId},
#'   \code{bucketId}, \code{bucketName} and \code{bucketType} will all be echoed
#'   back to the user.
#'
#' @examples
#' \dontrun{
#' # Make a bucket Private
#' b2UpdateBucket(bucketId = "aUniqueBucketId", bucketType = "allPrivate")
#' # Make a bucket Public
#' b2UpdateBucket(bucketId = "aUniqueBucketId", bucketType = "allPublic")
#' }
#'
#' @export

b2UpdateBucket <- function(bucketId, bucketType) {
  # Read Account Authorisation file
  accountAuthorization <- NULL
  accountAuthorization <- readRDS("accountAuthorization.rds")

  # Function options from input, make a dataframe
  accountId <- as.character(accountAuthorization$accountId)
  accountId <- as.data.frame(accountId, stringsAsFactors = FALSE)
  bucketId <- as.data.frame(bucketId, stringsAsFactors = FALSE)
  bucketType <- as.data.frame(bucketType, stringsAsFactors = FALSE)

  # Bind function option data frames together
  vars <- cbind(accountId, bucketId, bucketType)

  # API call
  b2Return <-
    httr::POST(
      paste(
        accountAuthorization$apiUrl,"/b2api/v1/b2_update_bucket", sep = ""
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
