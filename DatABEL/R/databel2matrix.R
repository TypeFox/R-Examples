#' converts 'databel' to matrix
#'
#' Converts a \code{\linkS4class{databel}} object to a regular R matrix.
#' This is the procedure used by the "as" converting to DatABEL objects,
#' in which case a temporary file name is created.
#'
#' @param from 'databel' matrix
#' @param rows which rows to include
#' @param cols which columns to include
#'
#' @return object of \code{\linkS4class{matrix}} class
#'
#' @author Stepan Yakovenko
#' @export
#'

databel2matrix <- function(from, rows, cols) {
    newi <- convert_intlogcha_index_to_int(rows,from,1)
    newj <- convert_intlogcha_index_to_int(cols,from,2)

    if (length(newi) * length(newj) > (2^31 -1))
      stop("The matrix has too many elements for R, see ?'Memory-limits'.")

    ret <- matrix(1.2345678912345, length(newi),length(newj));

    if(!.Call("assignDoubleMatrix", from@data, newi, newj, ret, as.integer(1)))
            stop("databel [<-: can't write variable.");

    dimnames(ret) <- dimnames(from)
        return(ret)
}
