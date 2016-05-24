#' @title
#' Prints an horizontal rule.
#'
#' @description
#' \code{hr}, given a list of symbols, composes with them an horizontal rule that fits the current width of the terminal window.
#' If nothing is passed to it, this function prints the default character (see \code{hrr.symbol} option) till the end of your current terminal window.
#' The same happens if symbol is an empty string or a string containing only whitespaces.
#'
#' @param       ... 		Symbols that compose the horizontal rule
#' @return      A list (invisibly).
#' @examples
#' hr()
#' hr('*')
#' hr('-', '#', '-')
#' hr('', '?', '   ')
#' @export
hr <- function(...) {
  symbols <- list(...)
  if (length(symbols) == 0) {
  	symbols <- list(getOption('hrr.symbol'))
  }
  symbols <- as.list(gsub("^\\s*$", getOption('hrr.symbol'), symbols))
  ncols <- ncols(set_option = FALSE)
  invisible(sapply(symbols, function(symbol) {
    repeat_count <- as.integer(ceiling(as.double(ncols) / nchar(symbol)))
    cat(substr(paste(rep_len(symbol, repeat_count), collapse = ''), 1, ncols), sep = '', '\n')
  }))
}
