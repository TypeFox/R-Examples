#' List B2 File Versions.
#'
#' \code{b2ListFileVersions} lists all versions of all of the files contained in
#' one bucket, within a user's account on the Backblaze B2 cloud storage
#' product.
#'
#' This function lists all versions of all of the files contained in one bucket,
#' within a user's account on the Backblaze B2 cloud storage product. Files will
#' be listed in alphabetical order by file name, and by reverse of date/time
#' uploaded for versions of files with the same name. Further details regarding
#' this API call are available here:
#'
#' \url{https://www.backblaze.com/b2/docs/b2_list_file_versions.html}
#'
#' \code{bucketId} is mandatory and must be user defined. \code{startFileName},
#' \code{startFileId} and \code{maxFileCount} are optional and may be defined by the
#' user if so desired.
#'
#' @param bucketId The unique identifier of the bucket containing the files to
#'   be listed. Bucket IDs may be obtained through the
#'   \code{b2ListBuckets}function in this package.
#' @param startFileName The name of the file from which the list will start. If
#'   there are no files with this name, the first version of the file with the
#'   first name after the given name will be the first in the list. This is an
#'   optional parameter. Not defining this parameter will result in the list
#'   starting from the newest file first. File names may be obtained through the
#'   \code{b2ListFileNames} function in this package.
#' @param startFileId The first file ID to return. This is an optional
#'   parameter. \code{startFileName} must also be provided if \code{startFileId}
#'   is specified.
#' @param maxFileCount An integer defining the maximum number of file names to
#'   return. This is an optional parameter and defaults to 100. The maximum
#'   acceptable value is 1000.
#' @return If successful a list will be returned containing \code{files},
#'   \code{nextFileName} and \code{nextFileId} for files within the specified
#'   bucket. If greater than the maximum number of specified files in
#'   \code{maxFileCount} exists, further file versions may be obtained beginning
#'   with \code{nextFileName}, as the \code{startFileName} parameter value. This
#'   may be done either alone or in combination with \code{nextFileId}, as
#'   \code{startFileId}. If successful a data frame will be returned, nested
#'   within \code{files}, containing \code{fileId}, \code{fileName},
#'   \code{action}, \code{size} and \code{uploadTimestamp} for all file versions
#'   within the specified bucket.
#'
#' @examples
#' \dontrun{
#' b2ListFileVersions(bucketId = "aUniqueBucketId",
#' startFileName = "yourFileName.png",
#' startFileId = "yourFileId"
#' maxFileCount = 500)
#' }
#'
#' @export

b2ListFileVersions <-
  function(bucketId, startFileName = "", startFileId = "", maxFileCount = 100) {
    # Read Account Authorisation file
    accountAuthorization <- NULL
    accountAuthorization <- readRDS("accountAuthorization.rds")

    # Function options from input, make a dataframe
    bucketId <- as.data.frame(bucketId, stringsAsFactors = FALSE)
    startFileName <-
      as.data.frame(startFileName, stringsAsFactors = FALSE)
    startFileId <-
      as.data.frame(startFileId, stringsAsFactors = FALSE)
    maxFileCount <-
      as.data.frame(maxFileCount, stringsAsFactors = FALSE)

    # Bind function option data frames together
    vars <- cbind(bucketId, startFileName, startFileId, maxFileCount)
    # Remove empty values
    vars <- vars[, colSums(vars != "") != 0]

    # API call
    b2Return <-
      httr::POST(
        paste(
          accountAuthorization$apiUrl,"/b2api/v1/b2_list_file_versions", sep = ""
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
