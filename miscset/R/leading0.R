#' @name leading0
#' @keywords zero
#' @author Sven E. Templer \email{sven.templer@@gmail.com}
#' @title Numeric to Character with Leading Zero(s)
#' @description 
#' Transform numeric values to character string prepending leading zero(s).
#' @param num Numeric vector (character also possible) to transform.
#' @param digits Numeric value of minimum length of output strings.
#' @return
#' Character vector with same length of strings of each value. 
#' Original "string" is prepended by zero(s). 
#' String length is at least \code{max(nchar(as.character(num)))}.
#' @examples
#' #
#' 
#' # use with paste to generate strings of equal size:
#' paste0("observation", leading0(1:10, 3))
#' 
#' #

#' @export leading0
leading0 <- function (num, digits = 2) {
  m <- nchar(num)
  digits <- max(digits, max(m))
  zeros <- paste(rep(0, digits), collapse="")
  paste0(substring(zeros, 0, digits - m), num)
}

# Many thanks to Bill Venables for the improved version.
