#' Upload File to B2.
#'
#' \code{b2Uploadfile} Uploads one file to the user's account B2.
#'
#' This function uploads one file to the user's account on the Backblaze B2
#' cloud storage product, returning a unique file ID. Files cannot be uploaded
#' to B2 unless \code{GetUploadUrl} has been executed first, specifying the
#' bucket in which to upload the file. Further details regarding this API call
#' are available here:
#'
#' \url{https://www.backblaze.com/b2/docs/b2_upload_file.html}
#'
#' \code{authToken}, \code{uploadUrl} and \code{fileName} are mandatory and must
#' be user defined.
#'
#' @param authToken A unique authorisation token permiting the upload of files
#'   to B2. An authorisation token is obtained from \code{b2UploadUrl}.
#' @param uploadUrl A unique URL permiting the upload of files to B2. An upload
#'   URL is obtained from \code{b2UploadUrl}.
#' @param fileName The name of the file to be uploaded. Files must be located in
#'   the current working directory. A complete list of supported file types may
#'   be found in the b2FileTypes.rds file located in this package's data
#'   directory, or at the following location:
#'   \url{https://www.backblaze.com/b2/docs/content-types.html}
#' @return If successful a list containing the \code{fileId}, \code{fileName},
#'   \code{accountId}, \code{bucketId}, \code{contentLength},
#'   \code{contentSha1}, \code{contentType} and \code{fileInfo} will all be
#'   echoed back to the user.
#'
#' @examples
#' \dontrun{
#' # Make a bucket Private
#' # Get Upload URL
#' uploadUrlReturn <- b2GetUploadUrl(bucketId = "aUniqueBucketId")
#' uploadUrl <- uploadUrlReturn$uploadUrl
#' authToken <- uploadUrlReturn$authorizationToken
#' # Upload file
#' b2UploadFile(authToken, uploadUrl, fileName = "yourFileName.png")
#' }
#'
#' @export

b2UploadFile <- function(authToken, uploadUrl, fileName) {
  # Read Account Authorisation file
  accountAuthorization <- NULL
  accountAuthorization <- readRDS("accountAuthorization.rds")

  # Function options from input, make a dataframe
  # File Name
  fileNameEncoded <- utils::URLencode(fileName, reserved = TRUE)
  # File Size
  fileSize <- file.size(fileName)
  # sha1 hash of file
  sha1Hash <- openssl::sha1(file(fileName))
  sha1Hash <- paste(sha1Hash, collapse = "")

  # Get file extension
  fileNameExtension <- tools::file_ext(fileName)
  # File Types List
  b2FileTypesLocation <- system.file("extdata", "b2FileTypes.rds", package = "backblazer")
  b2FileTypes <- readRDS(b2FileTypesLocation)
  b2ContentType <-
    b2FileTypes$contentType[grepl(fileNameExtension, b2FileTypes$fileExtension) == TRUE]

  # API call
  b2Return <-
    httr::POST(
      uploadUrl, body = httr::upload_file(fileName), httr::add_headers(
        'Authorization' = as.character(authToken), 'X-Bz-File-Name' = as.character(fileNameEncoded), 'Content-Type' = as.character(b2ContentType), 'Content-Length' = as.character(fileSize), 'X-Bz-Content-Sha1' = as.character(sha1Hash)
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
