#' Exports DatABEL object to a text file
#'
#' Exports DatABEL object to a text file
#'
#' @param databel DatABEL object
#' @param file output file name
#' @param NAString string to replace NA with
#' @param row.names export row names if TRUE
#' @param col.names export col names if TRUE
#' @param transpose whether the matrix should be transposed
#'
#' @author Stepan Yakovenko
#' @export
#'

databel2text <- function(databel, file, NAString = "NA",
                         row.names=TRUE, col.names=TRUE, transpose =
                         FALSE) {
    if (!is.character(file)) stop("databel save_as: file argument should be character")
    if (!.Call("saveAsText", databel@data, file,c(row.names,
                                                  col.names,transpose), NAString))
        stop("can not databel2text(): saveAsText failed")
    return(databel)
}
