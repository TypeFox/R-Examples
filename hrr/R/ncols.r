#' @title
#' Computes the maximum width of terminal window.
#'
#' @description
#' \code{\link{ncols}} detects the maximum number of columns of the current terminal window.
#' Note that columns are normally the same as characters except in CJK languages.
#'
#' @param     set_option    whether to automatically set the R terminal width to the system terminal maximum number of columns [default \code{TRUE}]
#' @return    The number of terminal window columns.
#' @export
ncols <- function(set_option = TRUE) {
  sysname <- Sys.info()['sysname']
  if ((ncol <- Sys.getenv('COLUMNS')) != '') {
    ncol <- as.integer(ncol)
    if (set_option) {
      options(width = ncol)
    }
    return(ncol)
  } else if ((sysname == 'Darwin' || sysname == 'Linux') && Sys.getenv('TERM') != '') {
    output <- tryCatch(system('tput cols', intern = TRUE), error = I)
    if (length(output) > 0) {
      ncol <- as.integer(sub('([0-9]+)', '\\1', output[1]))
      if (is.finite(ncol) && ncol > 0 && set_option) {
        options(width = ncol)
      }
    }
    rm(output)
    return(ncol)
  } else if (sysname == 'Windows') {
    return(getOption('width') - 1)
  }
  # fallback, should work also for Windows and various IDE (e.g. RStudio)
  return(getOption('width'))
}
