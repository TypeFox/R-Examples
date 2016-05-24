#' converts matrix to 'databel'
#'
#' Converts regular R matrix to \code{\linkS4class{databel}} object.
#' This is the procedure used by "as" converting to DatABEL objects,
#' in which case a temporary file name is created
#'
#' @param from R matrix
#' @param filename which FILEVECTOR BASE file name to use
#' @param cachesizeMb cache size to be used when accessing the object
#' @param type type of data to use for storage ("DOUBLE", "FLOAT", "INT",
#' "UNSIGNED_INT", "UNSIGNED_SHORT_INT", "SHORT_INT", "CHAR", "UNSIGNED_CHAR")
#' @param readonly whether to generate new 'databel' in read only mode
#'
#' @return object of class \code{\linkS4class{databel}}
#'
#' @author Yurii Aulchenko
#' @export
#'

matrix2databel <- function(from, filename, cachesizeMb = 64, type =
                           "DOUBLE", readonly = FALSE) {
    # a bit dirty
    to <- make_empty_fvf(as.character(filename),
                         nvariables=dim(from)[2],
                         nobservations=dim(from)[1],
                         type=type,
                         readonly=FALSE,
                         cachesizeMb=cachesizeMb)
    to[] <- from
    matrix_dimnames <- dimnames(from)
    if (!is.null(matrix_dimnames)) {
        set_dimnames(to) <- matrix_dimnames
    }
    setReadOnly(to) <- readonly
    return(to)
}
