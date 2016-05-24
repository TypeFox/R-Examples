#' @name strrev
#' @keywords inverse reverse
#' @author Sven E. Templer
#' @title Reverse Text Strings
#' @description 
#' Create a reverse version of strings.
#' @param x vector with strings. Is coerced to character.
#' @return
#' Returns a character vector with reversed strings.
#' @seealso
#' \link{rev}
#' @examples 
#' #
#' 
#' strrev(c("abc","asdf"))
#' 
#' #

#' @export strrev
strrev <- function (x) {
  
  xna <- is.na(x)
  x <- as.character(x)
  x <- strsplit(x, "")
  x <- lapply(x, function (y) {
    y <- paste(rev(y), collapse="")
    return(y)
  })
  x <- unlist(x)
  x[xna] <- NA
  return(x)

}
