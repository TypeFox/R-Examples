#' Download Dryad files
#'
#' @export
#' @param url Dryad URL for a dataset.
#' @param destfile Destination file. If not given, we assign a file name based
#' on URL provided.
#' @param ... Further args passed on to \code{\link{download.file}}
#' @return A path to the file
#' @details This function is a thin wrapper around download.file to get files
#' to your machine only. We don't attempt to read/parse them in to R.
#' @examples \dontrun{
#' url <- download_url('10255/dryad.1759')
#' dryad_fetch(url)
#' }
dryad_fetch <- function(url, destfile = NULL, ...) {
  if (is.null(destfile)) {
    num <- sub("/", "_", strextract(url, "10255/dryad\\.[0-9]+/[A-Za-z\\.]+"))
    destfile <- file.path(Sys.getenv("HOME"), num)
    if (!file.exists(destfile)) {
      dir.create(dirname(destfile), recursive = TRUE)
    }
  }
  download.file(url, destfile = destfile, ...)
}
