## methods for BerlinData generic functions ##

#' Gets the resources from a dataset
#' @export 
#' @method resources berlin_data_dataset
#' @param object a dataset with a list of resources
#' @param ... optional additional arguments
resources.berlin_data_dataset <- function(object, ...) {
  object$resources
}

#' @export
download.berlin_data_dataset <- function(x, ...) {
  message(paste("Downloading all resources for dataset: ", x$title))
  if (length(x$resources) == 0) {
    message(paste("  No resources located"))
    return()
  }
  y <- download(x$resources, ...)
  y
}

#' @export
getDatasetMetaData.berlin_data_dataset <- function(where, ...) where

## methods for base generic functions ##

#' @export
as.data.frame.berlin_data_dataset <- function(x, ...) {
  if (length(x$resources) == 0) {
    message(paste("No resources located for dataset: ", x$title))
    return()
  }
  y <- as.data.frame(x$resources, ...)
  y$dataset_title <- x$title
  y
}

#' @export
summary.berlin_data_dataset <- function(object, ...) {
  cat(object$title)
  cat("\n")
  cat(summary(object$resources))
}

#' @export
print.berlin_data_dataset <- function(x, ...) {
  summary(x, ...)
}

# roxygen2 doesn't recognize 'is.x' as an S3 method, requires manual documentation
#' @method berlin_data_dataset is
is.berlin_data_dataset <- function(x) inherits(x, "berlin_data_dataset")

