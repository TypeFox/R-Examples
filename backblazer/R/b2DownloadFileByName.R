#' Download B2 File by Name.
#'
#' \code{b2DownloadFileByName} downloads a file from the user's account on the
#' Backblaze B2 cloud storage product.
#'
#' This function downloads a file from the user's account on the Backblaze B2
#' cloud storage product using the file's name only. Files of the same name may
#' have multiple versions stored on B2. Therefore, only the most recent version
#' matching the specified filename will be downloaded. Further details regarding
#' this API call are available here:
#'
#' \url{https://www.backblaze.com/b2/docs/b2_download_file_by_name.html}
#'
#' \code{fileName} and \code{bucketName} are mandatory and must be user defined.
#' \code{overwrite} is optionally user defined and defaults to FALSE.
#'
#' @param fileName The name of the file to be downloaded. File names may be
#'   obtained through the \code{b2ListFileNames} function in this package.
#' @param bucketName The name of the bucket containing the requested file.
#'   Bucket names may be obtained through the \code{b2ListBuckets} function in
#'   this package.
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
#' b2DownloadFileByName(bucketName = "this-is-a-uniquely-named-bucket",
#' fileName = "yourFileName.txt", overwrite = TRUE)
#' }
#'
#' @export

b2DownloadFileByName <-
  function(bucketName, fileName, overwrite = FALSE) {
    # Read Account Authorisation file
    accountAuthorization <- NULL
    accountAuthorization <- readRDS("accountAuthorization.rds")

    # Function options from input, make a dataframe
    fileName <- as.data.frame(fileName, stringsAsFactors = FALSE)
    bucketName <- as.data.frame(bucketName, stringsAsFactors = FALSE)

    # Bind function option data frames together
    vars <- cbind(bucketName, fileName)

    # API call
    b2Return <-
      httr::GET(
        url = paste(
          accountAuthorization$downloadUrl,"/file/", bucketName, "/", fileName, sep =
            ""
        ), httr::add_headers(
          'Authorization' = as.character(accountAuthorization$authorizationToken)
        ), httr::write_disk("tmp", overwrite = overwrite)
      )

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
