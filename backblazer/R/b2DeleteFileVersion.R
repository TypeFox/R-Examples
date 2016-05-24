#' Delete B2 File Version.
#'
#' \code{b2DeleteFileVersion} deletes a version of a file in the user's account
#' on the Backblaze B2 cloud storage product.
#'
#' This function deletes a version of a file within the user's account on the
#' Backblaze B2 cloud storage product. Files of the same name may have multiple
#' versions stored on B2. If the deleted file version is the latest, and there
#' are older versions of the same file, the most recent of these will become the
#' latest version. The most recent file version is always the version downloaded
#' if requested by name. Further details regarding this API call are available
#' here:
#'
#' \url{https://www.backblaze.com/b2/docs/b2_delete_file_version.html}
#'
#' \code{fileName} and \code{fileId} are mandatory and must be user defined.
#'
#' @param fileName The name of the file to be deleted. File names may be
#'   obtained through the \code{b2ListFileNames} function in this package.
#' @param fileId The unique identifier of the file to be deleted. File IDs may
#'   be obtained through the \code{b2UploadFile}, \code{b2ListFileNames}, or
#'   \code{b2ListFileVersions} functions in this package.
#' @return If successful a list containing the \code{fileId} and \code{fileName}
#'   will be echoed back to the user.
#'
#' @examples
#' \dontrun{
#' b2DeleteFileVersion(fileName = "nameOfTheFileToDelete",
#'  fileId = "Unique_identifier_of_the_file_to_delete")
#' }
#'
#' @export

b2DeleteFileVersion <- function(fileName, fileId) {
  # Read Account Authorisation file
  accountAuthorization <- NULL
  accountAuthorization <- readRDS("accountAuthorization.rds")

  # Function options from input, make a dataframe
  fileName <- as.data.frame(fileName, stringsAsFactors = FALSE)
  fileId <- as.data.frame(fileId, stringsAsFactors = FALSE)

  # Bind function option data frames together
  vars <- cbind(fileName, fileId)

  # API call
  b2Return <-
    httr::POST(
      paste(
        accountAuthorization$apiUrl,"/b2api/v1/b2_delete_file_version", sep = ""
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
