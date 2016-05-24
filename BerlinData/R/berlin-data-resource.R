## methods for BerlinData generic functions ##

#' Downloads a resource
#' @param x a resource of type berlin_data_resource
#' @param message.on.fail logical: show message on download failure?
#' @param message.on.succeed logical: show message on sucessful download?
#' @param ... optional additional arguments to download function
#' @export
download.berlin_data_resource <- function(x, message.on.fail=TRUE, message.on.succeed=TRUE, ...) {
  resource_link <- x$url
  if(message.on.fail & substr(resource_link, 1, 5) == "https") {
    message(paste("Could not download resource:", resource_link))
    message("  R does not support the https:// URL scheme.")
    message("  Please try manual download.")
    return()
  }
  class(resource_link) <- x$format
  result <- download(resource_link, message.on.fail=message.on.fail, ...)
  if (message.on.succeed & length(result) > 0) message(paste("    Downloaded", x$format,":", resource_link))
  result
}

## methods for base generic functions ##

# roxygen2 doesn't recognize 'is.x' as an S3 method, requires manual documentation
#' @method berlin_data_resource is
is.berlin_data_resource <- function(x) inherits(x, "berlin_data_resource")

#' @export
summary.berlin_data_resource <- function(object, ...) {
  cat(paste0("Title: '", object$title,"'; Format: ", object$format, "\n"))
  cat(paste0("       URL: ", object$url))
}

#' @export
print.berlin_data_resource <- function(x, ...) {
  summary(x, ...)
}

#' @export
as.data.frame.berlin_data_resource <- function(x, ...) {
  y <- as.data.frame.list(x, ...)
  y$scheme <- factor(strsplit(x$url, '://')[[1]][1])
  y
}
