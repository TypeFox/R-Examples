#' Delete B2 Bucket.
#'
#' \code{b2DeleteBucket} deletes an existing bucket in the user's account on the
#' Backblaze B2 cloud storage product.
#'
#' This function deletes an existing bucket within the user's account on the
#' Backblaze B2 cloud storage product. Further details regarding this API call
#' are available here:
#'
#' \url{https://www.backblaze.com/b2/docs/b2_delete_bucket.html}
#'
#' \code{bucketId} is mandatory and must be user defined.
#'
#' @param bucketId The unique identifier of the bucket to be deleted. A
#'   list of all the user's bucket IDs may be found using the
#'   \code{b2_list_buckets} function in this package.
#' @return If successful a list containing the \code{accountId},
#'   \code{bucketId}, \code{bucketName} and \code{bucketType} of the deleted
#'   bucket will all be echoed back to the user.
#'
#' @examples
#' \dontrun{
#' b2DeleteBucket(bucketId = "aUniqueBucketId")
#' }
#'
#' @export

b2DeleteBucket <- function(bucketId) {
  # Read Account Authorisation file
  accountAuthorization <- NULL
  accountAuthorization <- readRDS("accountAuthorization.rds")

  # Function options from input, make a dataframe
  accountId <- as.character(accountAuthorization$accountId)
  bucketId <- as.data.frame(bucketId, stringsAsFactors = FALSE)

  # Bind function option data frames together
  vars <- cbind(accountId, bucketId)

  # API call
  b2Return <-
    httr::POST(
      paste(
        accountAuthorization$apiUrl,"/b2api/v1/b2_delete_bucket", sep = ""
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
