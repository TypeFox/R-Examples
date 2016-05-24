
#' Wait for a single keypress at the terminal
#'
#' It currently only works at Linux/Unix and OSX terminals,
#' see \code{\link{has_keypress_support}}.
#'
#' The following special keys are supported:
#' \itemize{
#'   \item Arrow keys: \sQuote{up}, \sQuote{down}, \sQuote{right},
#'     \sQuote{left}.
#'   \item Function keys: from \sQuote{f1} to \sQuote{f12}.
#'   \item Others: sQuote{home}, \sQuote{end},
#'     \sQuote{insert}, \sQuote{delete}, \sQuote{pageup},
#'     \sQuote{pagedown}.
#' }
#'
#' @return The key pressed, a character scalar.
#'
#' @family keypress
#' @useDynLib keypress
#' @export
#' @examples
#' \dontrun{
#' x <- keypress()
#' cat("You pressed key", x, "\n")
#' }

keypress <- function() {
  if (!has_keypress_support()) {
    stop("Your platform/terminal does not support keypress")
  }
  .Call("keypress", PACKAGE = "keypress")
}

#' Check if the current platform/terminal supports reading
#' single keys.
#'
#' @return Whether there is support for waiting for individual
#' keypressses.
#'
#' @family keypress
#' @export
#' @examples
#' has_keypress_support()

has_keypress_support <- function() {
  ## Supported if we have a terminal, and we are not in RStudio,
  ## not in R.app, not in Rgui, and not in Emacs.
  ## Yes, pretty limited.
  isatty(stdin()) &&
    Sys.getenv("RSTUDIO") != 1 &&
    Sys.getenv("R_GUI_APP_VERSION") == "" &&
    .Platform$GUI != "Rgui" &&
    ! identical(getOption("STERM"), "iESS") &&
    Sys.getenv("EMACS") != "t"
}
