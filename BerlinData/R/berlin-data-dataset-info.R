## methods for BerlinData generic functions ##

#' @export
getDatasetMetaData.berlin_data_dataset_info <- function(where, ...) {
  link <- where$link
  result <- parseMetaData(link, ...)
  result
}

#' @export
download.berlin_data_dataset_info <- function(x, ...) {
  message(paste(" Cannot download dataset:", x$title))
  message("  Please call getDatasetMetaData to get a list of available resources for download.")
}

## methods for base generic functions ##

# roxygen2 doesn't recognize 'is.x' as an S3 method, requires manual documentation
#' @method berlin_data_dataset_info is
is.berlin_data_dataset_info <- function(x) inherits(x, "berlin_data_dataset_info")

#' @export
as.data.frame.berlin_data_dataset_info <- function(x, ...) {
  y <- as.data.frame.list(x, ...)
  y
}
