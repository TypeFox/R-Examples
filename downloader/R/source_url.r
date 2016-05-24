#' Download an R file from a URL and source it
#'
#' This will download a file and source it. Because it uses the
#' \code{\link{download}()} function, it can handle https URLs.
#'
#' By default, \code{source_url()} checks the SHA-1 hash of the file. If it
#' differs from the expected value, it will throw an error. The default
#' expectation is that a hash is provided; if not, \code{source_url()} will
#' prompt the user, asking if they are sure they want to continue, unless
#' \code{prompt=FALSE} is used. In other words, if you use \code{prompt=FALSE},
#' it will run the remote code without checking the hash, and without asking
#' the user.
#'
#' The purpose of checking the hash is to ensure that the file has not changed.
#' If a \code{source_url} command with a hash is posted in a public forum, then
#' others who source the URL (with the hash) are guaranteed to run the same
#' code every time. This means that the author doesn't need to worry about the
#' security of the server hosting the file. It also means that the users don't
#' have to worry about the file being replaced with a damaged or
#' maliciously-modified version.
#'
#' To find the hash of a local file, use \code{\link{digest}()}. For a simple
#' way to find the hash of a remote file, use \code{\link{sha_url}()}.
#'
#' @param url The URL to download.
#' @param sha A SHA-1 hash of the file at the URL.
#' @param prompt Prompt the user if no value for \code{sha} is provided.
#' @param quiet If \code{FALSE} (the default), print out status messages about
#'   checking SHA.
#' @param ... Other arguments that are passed to \code{\link{source}()}.
#'
#' @seealso \code{\link{source}()} for more information on the arguments
#'   that can be used with this function. The \code{\link{sha_url}()} function
#'   can be used to find the SHA-1 hash of a remote file.
#'
#' @export
#' @examples
#' \dontrun{
#' # Source the a sample file
#'
#' # This is a very long URL; break it up so it can be seen more easily in the
#' # examples.
#' test_url <- paste0("https://gist.github.com/wch/dae7c106ee99fe1fdfe7",
#'                    "/raw/db0c9bfe0de85d15c60b0b9bf22403c0f5e1fb15/test.r")
#' downloader::source_url(test_url,
#'                        sha = "9b8ff5213e32a871d6cb95cce0bed35c53307f61")
#'
#' # Find the hash of a file
#' downloader::sha_url(test_url)
#' }
#'
#'
#' @importFrom digest digest
source_url <- function(url, sha = NULL, ..., prompt = TRUE, quiet = FALSE) {

    if (prompt && (is.null(sha) || sha == '')) {
      resp <- readline(prompt = paste(sep = '',
        ' No SHA-1 hash specified for the file. The hash is needed to ensure that\n',
        ' the file at the URL has not changed. See ?source_url for information on\n',
        ' why this is useful. Are sure you want to continue? [y/n] '))

      sha <- NULL  # Set to NULL for simpler check later on

      if (tolower(resp) != "y") {
        message("Quitting")
        return(invisible())
      }
    }

    temp_file <- tempfile()
    download(url, temp_file)
    on.exit(unlink(temp_file))

    if (!is.null(sha)) {
      url_sha <- digest(file = temp_file, algo = 'sha1')

      if (url_sha == sha) {
        if (!quiet) {
          message('Hash ', url_sha, ' matches expected value.')
        }
      } else {
        stop('Hash ', url_sha, ' does not match expected value!')
      }

    } else {
      if (!quiet) {
        message('Not checking SHA-1 of downloaded file.')
      }
    }

    source(temp_file, ...)
}
