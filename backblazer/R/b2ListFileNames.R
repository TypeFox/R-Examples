#' List B2 File Names.
#'
#' \code{b2ListFileNames} lists the names of all files in a bucket, starting at
#' a given name.
#'
#' This function lists the names of all files in a bucket, starting at a given
#' name, within a user's account on the Backblaze B2 cloud storage product.
#' There may be many file versions for the same name, but this function will
#' return each name only once. If you want all of the versions, use
#' \code{b2ListFileVersions} instead. Further details regarding this API call
#' are available here:
#'
#' \url{https://www.backblaze.com/b2/docs/b2_list_file_names.html}
#'
#' \code{bucketId} is mandatory and must be user defined. \code{startFileName}
#' and \code{maxFileCount} are optional and may be defined by the user if so
#' desired.
#'
#' @param bucketId The unique identifier of the bucket containing the files to
#'   be listed. Bucket IDs may be obtained through the
#'   \code{b2ListBuckets}function in this package.
#' @param startFileName The name of the file from which the list will start.
#'   This is an optional parameter. Not defining this parameter will result in
#'   the list starting from the newest file first. File names may be obtained
#'   through the \code{b2ListFileNames} function in this package.
#' @param maxFileCount An integer defining the maximum number of file names to
#'   return. This is an optional parameter and defaults to 100. The maximum
#'   acceptable value is 1000.
#' @return If successful a list will be returned containing \code{files} and
#'   \code{nextFileName}, for files within the specified bucket. If greater than
#'   the maximum number of specified files in \code{maxFileCount} exists,
#'   further file names may be obtained beginning with \code{nextFileName}, as
#'   the \code{startFileName} parameter value. If successful a data frame will
#'   be returned, nested within \code{files}, containing \code{fileId},
#'   \code{fileName}, \code{action}, \code{size} and \code{uploadTimestamp} for
#'   all files within the specified bucket.
#'
#' @examples
#' \dontrun{
#' b2ListFileNames(bucketId = "aUniqueBucketId",
#' startFileName = "yourFileName.png",
#' maxFileCount = 500)
#' }
#'
#' @export

b2ListFileNames <-
  function(bucketId, startFileName = "", maxFileCount = 100) {
    # Read Account Authorisation file
    accountAuthorization <- NULL
    accountAuthorization <- readRDS("accountAuthorization.rds")

    # Function options from input, make a dataframe
    bucketId <- as.data.frame(bucketId, stringsAsFactors = FALSE)
    startFileName <-
      as.data.frame(startFileName, stringsAsFactors = FALSE)
    maxFileCount <-
      as.data.frame(maxFileCount, stringsAsFactors = FALSE)

    # Bind function option data frames together
    vars <- cbind(bucketId, startFileName, maxFileCount)
    # Remove empty values
    vars <- vars[, colSums(vars != "") != 0]

    # API call
    b2Return <-
      httr::POST(
        paste(
          accountAuthorization$apiUrl,"/b2api/v1/b2_list_file_names", sep = ""
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
