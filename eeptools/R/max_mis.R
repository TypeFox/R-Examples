##' A function to safely take the maximum of a vector that could include only NAs.
##' @description When computing the maximum on arbitrary subsets of data, some 
##' of which may only have missing values, it may be necessary to take the maximum 
##' of a vector of NAs. This replaces the behavior that returns Inf or-Inf and 
##' replaces it with simply returning an NA. 
##' @param x A vector of data that a maximum can be taken of.
##' @return A vector with the maximum value or with an NA of the proper type
##' @details This function only returns valid results for vectors with a mix of 
##' NA and numeric values.
##' @seealso See also \code{\link{max}} which this function wraps.
##' @author Jared E. Knowles
##' @export
##' @examples
##' max(c(7,NA,3,2,0),na.rm=TRUE)
##' max_mis(c(7,NA,3,2,0))
##' max(c(NA,NA,NA,NA),na.rm=TRUE)
##' max_mis(c(NA,NA,NA,NA))
##' 
max_mis <- function(x){
  varclass <- class(x)
  if(varclass %in% c("ordered", "factor", "character")){
    stop("Vector must be of class integer or real to take maximum values")
  }
  suppressWarnings(x <- max(x, na.rm=TRUE))
  if(varclass == "integer"){
    ifelse(!is.finite(x), NA_integer_, x)
  } else if(varclass == "numeric") {
    ifelse(!is.finite(x), NA_real_, x)
  } else if(varclass %in% c("NULL", "logical")){
    NA_real_
  }
}


#' Find the nth maximum value
#'
#' @note If \code{n} is smaller/larger than \code{0}/\code{length(unique(x))}
#' the error \sQuote{index outside bounds} is thrown.
#'
#' @param x a vector of numeric values
#' @param n which max to return
#'
#' @return the value of the nth most maximum value in a vector
#' @export
#'
#' @examples
#' x <- c(1:20, 20:1)
#' nth_max(x, n = 1) #20
#' nth_max(x, n = 2) #19
nth_max <- function(x, n = 1){
  x <- unique(x)
  n <- length(x) - n + 1
  sort(x, partial=n)[n]
}
