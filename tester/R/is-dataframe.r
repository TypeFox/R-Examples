#' @title Is data frame
#' 
#' @description
#' \code{is_dataframe} tests if an object is a data frame \cr
#' \code{is_numeric_dataframe} tests if an object is a numeric data frame \cr
#' \code{is_string_dataframe} tests if an object is a string data frame \cr
#' \code{is_factor_dataframe} tests if an object is a factor data frame \cr
#' \code{is_not_dataframe} tests if an object is not a data frame 
#' 
#' @param x an R object
#' @name is_dataframe
#' @aliases is_dataframe is_numeric_dataframe is_string_dataframe
#' is_factor_dataframe is_not_dataframe
#' @export is_dataframe is_numeric_dataframe is_string_dataframe
#' is_factor_dataframe is_not_dataframe
#' @examples
#' is_dataframe(iris) # TRUE
#' is_dataframe(1:10) # FALSE
#' 
#' is_numeric_dataframe(iris) # FALSE
#' is_numeric_dataframe(iris[,1:4]) # TRUE
#' 
#' DF = matrix(letters[1:24], 6, 4)
#' DF1 = data.frame(DF)
#' DF2 = data.frame(DF, stringsAsFactors=FALSE)
#' 
#' is_string_dataframe(DF1) # FALSE
#' is_string_dataframe(DF2) # TRUE
#' 
#' is_factor_dataframe(DF1) # TRUE
#' is_factor_dataframe(DF2) # FALSE
NULL

is_dataframe <- function(x) {
  if (is.data.frame(x)) TRUE else FALSE
}

is_numeric_dataframe <- function(x) {
  if (is.data.frame(x)) {
    numerics = unlist(lapply(x, is.numeric))
    if (sum(numerics) == dim(x)[2L]) TRUE else FALSE
  } else FALSE
}

is_string_dataframe <- function(x) {
  if (is.data.frame(x)) {
    characters = unlist(lapply(x, is.character))
    if (sum(characters) == dim(x)[2L]) TRUE else FALSE
  } else FALSE
}

is_factor_dataframe <- function(x) {
  if (is.data.frame(x)) {
    factors = unlist(lapply(x, is.factor))
    if (sum(factors) == dim(x)[2L]) TRUE else FALSE
  } else FALSE
}

is_not_dataframe <- function(x) {
  !is_dataframe(x)
}
