#' Download a file from a URL and find a SHA-1 hash of it
#'
#' This will download a file and find a SHA-1 hash of it, using
#' \code{\link{digest}()}. The primary purpose of this function is to provide
#' an easy way to find the value of \code{sha} which can be passed to
#' \code{\link{source_url}()}.
#'
#' @param url The URL of the file to find a hash of.
#' @param cmd If \code{TRUE} (the default), print out a command for sourcing the
#'   URL with \code{\link{source_url}()}, including the hash.
#'
#' @export
#' @examples
#' \dontrun{
#' # Get the SHA hash of a file. It will print the text below and return
#' # the hash as a string. This is a very long URL; break it up so it can be
#' # seen more easily in the examples.
#' test_url <- paste0("https://gist.github.com/wch/dae7c106ee99fe1fdfe7",
#'                    "/raw/db0c9bfe0de85d15c60b0b9bf22403c0f5e1fb15/test.r")
#' sha_url(test_url)
#' # Command for sourcing the URL:
#' #  downloader::source_url("https://gist.github.com/wch/dae7c106ee99fe1fdfe7
#' #  /raw/db0c9bfe0de85d15c60b0b9bf22403c0f5e1fb15/test.r",
#' #    sha="9b8ff5213e32a871d6cb95cce0bed35c53307f61")
#' # [1] "9b8ff5213e32a871d6cb95cce0bed35c53307f61"
#' }
#'
#'
#' @importFrom digest digest
sha_url <- function(url, cmd = TRUE) {
  temp_file <- tempfile()
  download(url, temp_file)
  on.exit(unlink(temp_file))

  sha <- digest(file = temp_file, algo = 'sha1')

  if (cmd) {
    message('Command for sourcing the URL:\n',
      '  downloader::source_url("', url, '", sha="', sha, '")')
  }

  sha
}
