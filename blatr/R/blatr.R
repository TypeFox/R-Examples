#' 'Blat' is a feature-rich command line email tool for Windows. The
#' \code{blatr} package is a wrapper for using 'Blat' from within R. For
#' information about 'Blat' and its options, see \code{http://www.blat.net}, for
#' the 'Blat' license, see \code{http://www.blat.net/?docs/license.txt}.
#'
#' To use the \code{blatr} package, 'Blat' needs to be installed. This can be
#' done manually, by downloading the executables from the website
#' \code{www.blat.net} and place them in the \code{blatr} installation
#' directory, i.e. in \code{/path/to/R/library/blatr/}. You can also use
#' \code{blatr:::install_blat} to let \code{blatr} do it for you. You can always
#' check that the needed files are placed correctly using
#' \code{blatr:::check_install}. Note that these functions are not exported, as
#' they are generally not needed when \code{blatr} is properly set up. The files
#' installed by \code{blatr:::install_blat} are hosted on the blatr development
#' page at \url{http://GitHub.com/smbache/blatr}. To get started, see the
#' documentation for \code{\link{blat}}.
#'
#' @docType package
#' @name blatr
#' @title Send Emails Using 'Blat' for Windows
#' @author Stefan Milton Bache
NULL

#' Send Emails Using 'Blat' for Windows
#'
#' This is a wrapper which calls 'Blat' to send emails. For documentation on the
#' options, see \url{http://www.blat.net/syntax/syntax.html}. You should not use
#' dashes as part of the argument names (as with 'Blat'), this is done for you:
#' simply use \code{name = value} pairs as arguments. For some basic examples of
#' the most commonly used arguments, see below.
#'
#' 'Blat' can use a file as body of the email, in which case this is ordinarily
#' the first argument in the 'Blat' command. However, \code{blatr} uses the
#' named argument \code{filename} to specify this. Unnamed arguments do not work
#' with blatr.
#'
#' @param ... arguments to blat
#' @rdname blat
#' @export
#' @examples
#' \dontrun{
#'    # With attachment
#'    blat(f      = "Your Name <your@@email.com>",
#'         to     = "your-recipients@@email.com",
#'         s      = "The subject",
#'         server = "server.address",
#'         attach = "C:/path/to/attachment.txt",
#'         body   = "The text you wish to send.")
#'
#'    # With file as body
#'    blat(f        = "Your Name <your@@email.com>",
#'         to       = "your-recipients@@email.com",
#'         s        = "The subject",
#'         server   = "server.address",
#'         attach   = "C:/path/to/attachment.txt",
#'         filename = "C:/path/to/file/with/body.txt")
#'
#'    # With username and password required.
#'    blat(f        = "Your Name <your@@email.com>",
#'         to       = "your-recipients@@email.com",
#'         s        = "The subject",
#'         server   = "server.address",
#'         attach   = "C:/path/to/attachment.txt",
#'         filename = "C:/path/to/file/with/body.txt",
#'         u        = "username",
#'         pw       = "password")
#' }
blat <- function(...)
{
  # The blat files should exist! They are not bundled with the package, but can be
  # installed with blatr:::install_blat
  # Bundling binaries is not allowed by CRAN.
  if (!check_install())
    stop("Blat files not installed. Use blatr:::install_blat before first use.")

  # unlist the arguments.
  args <- list(...)

  # The filename argument is unnamed in blat, so deal with this explicitely.
  if ("filename" %in% names(args)) {
    # Get it and remove it from the list. It will always be added below, but
    # may be the empty string.
    filename <- shQuote(args[["filename"]])
    args <- args[names(args) != "filename"]
  } else {
    filename <- ""
    if (!"body" %in% names(args) || args[["body"]] == "")
      args[["body"]] <- "\n"
  }

  # Construct the blat command.
  cmd <-
    paste(
      lapply(seq_along(args), function(i) {
        arg <- shQuote(args[[i]])
        paste0("-", names(args)[i], " ", arg)
      }),
      collapse = " ")

  # prepend the name of the executable
  cmd <- paste(shQuote(blat_files()[1L]), filename, cmd)

  # Execute the command.
  system(cmd)
}
