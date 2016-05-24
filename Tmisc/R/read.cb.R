#' Read from the clipboard
#' 
#' Read tabular data from the clipboard.
#' 
#' @author Stephen Turner
#' @keywords keywords
#' 
#' @param header A logical value indicating whether the file contains the names of the variables as its first line. Overrides the default \code{header=FALSE} option in \code{read.table()}.
#' @param ... Further arguments to be passed to \code{read.table}
#' 
#' @return A data.frame
#' @import utils
#' 
#' @examples
#' \dontrun{
#' # To read CSV data with a header from the clipboard:
#' read.cb(header=TRUE, sep=',')
#' }
#' 
#' @export
read.cb <- function(header=TRUE, ...) {
    os <- Sys.info()[1]
    if (os=="Darwin")       content <- read.table(pipe("pbpaste"),  header=header, ...)
    else if (os=="Windows") content <- read.table(file="clipboard", header=header, ...)
    else stop("Sorry, only works on Windows and Mac.")
    if ("package:dplyr" %in% search()) content <- tbl_df(content)
    return(content)
}