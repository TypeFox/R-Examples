#' Create tables in dokuwiki format
#' 
#' Prints the supplied data frame or matrix using Dokuwiki's table syntax, optionally copying the data to the clipboard (Mac OS X only).
#' 
#' @author Stephen Turner
#' 
#' @param x A data.frame.
#' @param headersep The separator used between entries in the header row.
#' @param sep The separator used between entries in all other rows.
#' @param clip Whether or not to write the returned table to the clipboard (currently only supported on Mac OS X). 
#' @param ... Further arguments passed to \code{write.table}.
#' 
#' @examples
#' dokuwiki(head(iris), clip=FALSE)
#' dokuwiki(head(mtcars), clip=FALSE, row.names=TRUE)
#' 
#' @export
dokuwiki <- function(x, headersep="^", sep="|", clip=TRUE, ...) {
    .dots <- list(...)
    .dots$x <- x
    .dots$sep <- sep
    .dots$col.names <- FALSE
    .dots$row.names <- ifelse(is.null(.dots$row.names), FALSE, .dots$row.names)
    .dots$quote <- ifelse(is.null(.dots$quote), FALSE, .dots$quote)
    .dots$na <- ifelse(is.null(.dots$na), "", .dots$na)
    # Header row. If printing row.names, add an extra header separator column
    if (.dots$row.names) {
        row1 <- paste0(headersep, " ", headersep, paste(colnames(x), collapse=headersep), headersep, "\n")
    } else {
        row1 <- paste0(headersep, paste(colnames(x), collapse=headersep), headersep, "\n")
    }
    # All other rows
    otherrows <- paste0(sep, utils::capture.output(do.call(utils::write.table, .dots)), sep, collapse = "\n")
    allrows <- paste0(row1, otherrows, collapse="\n")
    if (clip) {
        if (Sys.info()["sysname"]=="Darwin") {
            con <- pipe("pbcopy")
            writeChar(allrows, con=con, eos=NULL)
            close(con)
            message("Copied to clipboard:\n")
        } else {
            warning("Writing to clipboard is supported on Mac OS X only.")
        }
    }
    cat(allrows)
}
