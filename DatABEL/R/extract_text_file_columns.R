#'
#' extracts columns from text file
#'
#' Extracts a column from text file to a matrix.
#' If in a particular file line the number of columns is less
#' then a column specified, returns last column!
#'
#' @param file file name
#' @param whichcols which columns to extract
#'
#' @return matrix of strings with values from that columns
#' @export
#'

extract_text_file_columns <- function(file, whichcols)
{
    if (!file.exists(file)) stop("file dose not exist")
    if (any(as.integer(whichcols) <= 0)) stop("whichcols must be positive integer")
    if (any(as.integer(whichcols) > 1000)) {
        warning("some whichcols >1000, are you sure?")
    }
    str <- .Call("extract_text_file_column_cpp",
                 as.character(file),
                 as.integer(whichcols-1))

    if (length(whichcols) > 1)
        {
            str <- matrix(str, nrow=length(whichcols))
            str <- t(str)
        }
    return(str)
}
