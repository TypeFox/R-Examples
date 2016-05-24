#' Hide B2 File.
#'
#' \code{b2HideFile} hides a file in the user's account on the Backblaze B2
#' cloud storage product.
#'
#' This function hides a file in the user's account on the Backblaze B2 cloud
#' storage product, so that it cannot be downloaded by name. In order to further
#' understand the concept of B2 file versions, see the speficic Backblaze
#' documenation. Further details regarding this API call are available here:
#'
#' \url{https://www.backblaze.com/b2/docs/b2_hide_file.html}
#' \url{https://www.backblaze.com/b2/docs/file_versions.html}
#'
#' \code{bucketId} \code{fileName} are mandatory and must be user defined.
#'
#' @param bucketId The unique identifier of the bucket containing the file to be
#'   hidden. Bucket IDs may be obtained through the \code{b2ListBuckets}
#'   function in this package.
#' @param fileName The name of the file to be hidden. File names may be obtained
#'   through the \code{b2ListFileNames} function in this package.
#' @return If successful a list will be returned containing \code{fileId},
#'   \code{fileName}, \code{action}, \code{size} and \code{uploadTimestamp}.
#'
#' @examples
#' \dontrun{
#' b2HideFile(bucketId = "aUniqueBucketId", fileName = "yourFileName.txt")
#' }
#'
#' @export

b2HideFile <- function(bucketId, fileName) {
  # Read Account Authorisation file
  accountAuthorization <- NULL
  accountAuthorization <- readRDS("accountAuthorization.rds")

  # Function options from input, make a dataframe
  fileName <- as.data.frame(fileName, stringsAsFactors = FALSE)
  bucketId <- as.data.frame(bucketId, stringsAsFactors = FALSE)

  # Bind function option data frames together
  vars <- cbind(bucketId, fileName)

  # API call
  b2Return <-
    httr::POST(
      paste(
        accountAuthorization$apiUrl,"/b2api/v1/b2_hide_file", sep = ""
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
