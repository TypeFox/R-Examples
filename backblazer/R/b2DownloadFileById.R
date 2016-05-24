#' Download B2 File by ID.
#'
#' \code{b2DownloadFileById} downloads a file from the user's account on the
#' Backblaze B2 cloud storage product.
#'
#' This function downloads a file from the user's account on the Backblaze B2
#' cloud storage product using the file's unique ID only. Files of the same name
#' may have multiple versions stored on B2. Therefore, every file version will
#' have a unique ID, which can be used for downloading that specific version.
#' Further details regarding this API call are available here:
#'
#' \url{https://www.backblaze.com/b2/docs/b2_download_file_by_id.html}
#'
#' \code{fileId} is mandatory and must be user defined. \code{overwrite} is
#' optionally user defined and defaults to FALSE.
#'
#' @param fileId The unique identifier of the file to be downloaded. File IDs
#'   may be obtained through the \code{b2ListFiles}, \code{b2ListFileVersions}
#'   and \code{b2UploadFile} functions in this package.
#' @param overwrite Binary TRUE or FALSE decision to overwrite any files in the
#'   current working directory, whose names match the downloaded file name.
#' @return If successful the response headers include the Content-Type that was
#'   specified when the file was uploaded. They also include the X-Bz-FileName
#'   and X-Bz-Content-Sha1 headers. The X-Bz-FileName uses percent-encoding, as
#'   if it were a URL parameter. If successful, the file will be downloaded to
#'   the current working directory.
#'
#' @examples
#' \dontrun{
#' b2DownloadFileById(fileId = "Unique_identifier_of_the_file_to_download",
#' overwrite = TRUE)
#' }
#'
#' @export

b2DownloadFileById <- function(fileId, overwrite = FALSE) {
  # Read Account Authorisation file
  accountAuthorization <- NULL
  accountAuthorization <- readRDS("accountAuthorization.rds")

  # Function options from input, make a dataframe
  fileId <- as.data.frame(fileId, stringsAsFactors = FALSE)

  # API call
  b2Return <-
    httr::POST(
      paste(
        accountAuthorization$downloadUrl,"/b2api/v1/b2_download_file_by_id", sep =
          ""
      ), body = jsonlite::toJSON(jsonlite::unbox(fileId), pretty = TRUE), httr::add_headers(
        'Authorization' = as.character(accountAuthorization$authorizationToken)
      ), httr::write_disk("tmp", overwrite = overwrite)
    )

  # Alternative GET call
  # b2Return <- httr::GET(url = paste(accountAuthorization$downloadUrl,"/b2api/v1/b2_download_file_by_id?fileId=", fileId, sep=""), add_headers('Authorization' = as.character(accountAuthorization$authorizationToken)), write_disk("tmp", overwrite))

  # Check for bad authorisation and sent message
  if (httr::status_code(b2Return) != "200") {
    badReturn <- jsonlite::fromJSON(httr::content(b2Return,type = "text"))
    stop(
      "\nSomething went wrong. Please check the function options to ensure valid values. \n",
      "\nStatus Code: ", badReturn$code, "\nMessage: ", badReturn$message
    )

  } else {
    # Rename tmp
    if (file.exists(b2Return$headers$'x-bz-file-name') &
        !isTRUE(overwrite)) {
      print("Unable to write to disk. File(s) exist and overwrite is set to FALSE")
    }

    else {
      renameResult <-
        file.rename(from = "tmp", to = b2Return$headers$'x-bz-file-name')
      # Output message
      print("File(s) downloaded successfully and saved to disk.")
    }
  }
}
