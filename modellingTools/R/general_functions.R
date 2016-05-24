# General helper functions

# Contents:
# - Replacement of the + operator to do automatic string concatenation
# - function for returning a vector from a column of a data frame

#==============================================================================#

# Replacement of the + operator

`+` <- function(a,b) {
  if (is.character(a) || is.character(b)) {
    return (stringr::str_c(a,b))
  }
  else {
    .Primitive('+')(a,b)
  }

}

#' Get the contents of a single column of a tbl
#'
#' Function to return a vector with the contents of one column of a tbl_df
#' The type of the column will be preserved if numeric or factor; all else will
#' be conveted to a vector of character strings
#'
#' @param  d a tbl
#' @param  i a character string containing the name of the column to be
#'         returned OR
#'         a single number representing the position of the column to be
#'         returned
#' @return  A vector containing the values from the specified column.
#'          The type will be numeric, factor, or character, depending on the
#'          type of the original column
#' @examples
#' x <- dplyr::tbl_df(iris)
#' column_vector(x,2)
#' column_vector(x,"Species")
#' @export

column_vector <- function(d,i) {
  if (length(i) > 1) warning("Requested column index is a vector containing" +
                              " multiple items. This may result in" +
                              " unexpected behaviour")
  col <- unlist(d[,i])
  if (is.numeric(col)) {
    return(as.numeric(col))
  } else if (is.factor(col)) {
    return(unname(as.factor(col)))
  } else {
    return(as.character(col))
  }
}
