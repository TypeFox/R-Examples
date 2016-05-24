#' Get B2 Upload URL.
#'
#' \code{b2GetUploadUrl} returns the URL required in order to upload files to
#' the user's account on the Backblaze B2 cloud storage product.
#'
#' This function returns the URL required in order to upload files the user's
#' account on the Backblaze B2 cloud storage product. An uploadUrl and upload
#' authorizationToken are also returned. These are valid for 24 hours or until
#' the endpoint rejects an upload, Further details regarding this API call are
#' available here:
#'
#' \url{https://www.backblaze.com/b2/docs/b2_get_upload_url.html}
#'
#' \code{bucketId} is mandatory and must be user defined.
#'
#' @param bucketId The unique identifier of the bucket where files are to be
#'   uploaded. Bucket IDs may be obtained through the
#'   \code{b2ListBuckets}function in this package.
#' @return If successful a list will be returned containing \code{bucketId},
#'   \code{uploadURL} and \code{authorizationToken}.
#'
#' @examples
#' \dontrun{
#' uploadUrlReturn <- b2GetUploadUrl(bucketId = "aUniqueBucketId")
#' uploadUrl <- uploadUrlReturn$uploadUrl
#' authToken <- uploadUrlReturn$authorizationToken
#' }
#'
#' @export

b2GetUploadUrl <- function(bucketId) {
  # Read Account Authorisation file
  accountAuthorization <- NULL
  accountAuthorization <- readRDS("accountAuthorization.rds")

  # Function options from input, make a dataframe
  bucketId <- as.data.frame(bucketId, stringsAsFactors = FALSE)

  # API call
  b2Return <-
    httr::POST(
      paste(
        accountAuthorization$apiUrl,"/b2api/v1/b2_get_upload_url", sep = ""
      ), body = jsonlite::toJSON(jsonlite::unbox(bucketId), pretty = TRUE), httr::add_headers(
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
