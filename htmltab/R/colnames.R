#' Produces correct column names
#'
#' @param df the generated body
#' @param header.names the header column names vector
#' @param colNames either a self-specificed character vector for the column names or a function used on header.names
#' @param xpath generated header and body xpath
#' @return a character vector of header column names
#' @noRd
make_colnames <- function(df, header.names = NULL, colNames = NULL, header.xpath) {

  if(length(header.names) == 0 && is.null(colNames)){
    warning("No header generated. Try passing information to header or colNames", call. = FALSE)
    return(df)
    }

  if(length(header.names) > 0 && is.null(colNames)) {
    if(ncol(df) != length(header.names)){
      warning("Header dimension doesn't match body dimension", call. = FALSE)
      colnames(df) <- vector()
    }
    colnames(df) <- header.names
    return(df)
  }

  if(is.character(colNames)) {
    header.names <- colNames
    if(ncol(df) != length(header.names)){
      warning("Header dimension doesn't match body dimension", call. = FALSE)
      colnames(df) <- vector()
    }
    colnames(df) <- header.names
    return(df)
  }

  if(class(colNames) == "function") {
    header.names <- colNames(header.names)
    if(ncol(df) != length(header.names)){
      warning("Header dimension doesn't match body dimension", call. = FALSE)
      colnames(df) <- vector()
    }
    colnames(df) <- header.names
  }

}
