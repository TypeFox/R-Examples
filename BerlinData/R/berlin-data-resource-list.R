## methods for BerlinData generic functions ##

#' @export
download.berlin_data_resource_list <- function(x, ...) {
  message(paste0(" Downloading ", length(x), " resource", ifelse(length(x)>1, "s", "")))
  y <- lapply(x, function(resource) tryCatch(download(resource, ...),
                                             error = function(e) {
                                               message(paste('  Could not download resource:', resource$title))
                                               message(paste('    ',e,'\n'))
                                               return()
                                             }))
  successful.downloads <- Filter(function(r)!is.null(r), y)
  message(paste(" Downloaded", length(successful.downloads), "of", length(x), "resources."))
  successful.downloads
}

## methods for base generic functions ##

#' @export
summary.berlin_data_resource_list <- function(object, ...) {
  cat(paste0(length(object), " resource", ifelse(length(object) > 1, "s", "")))
  cat("\n")
  for (i in 1:length(object)) {
    cat(paste0(i,": "))
    cat(summary(object[[i]]))
    cat("\n")
  }
}

#' @export
print.berlin_data_resource_list <- function(x, ...) {
  summary(x, ...)
}

#' @export
as.data.frame.berlin_data_resource_list <- function(x, ...) {
  y <- lapply(x, as.data.frame, ...)
  resource_list_lengths <- unlist(sapply(y, nrow))
  y <- do.call(rbind, y)
  y$position_in_resource_list <- rep(1:length(resource_list_lengths), times=resource_list_lengths)
  y
}

# roxygen2 doesn't recognize 'is.x' as an S3 method, requires manual documentation
#' @method berlin_data_resource_list is
is.berlin_data_resource_list <- function(x) inherits(x, "berlin_data_resource_list")


# subsets inherit class 
#' @export
`[.berlin_data_resource_list` <- function(x, ...) {
  x <- NextMethod("[")
  class(x) = "berlin_data_resource_list"
  x
}
