#' @name strextr
#' @keywords extract string
#' @author Sven E. Templer
#' @title Extract a Substring
#' @description 
#' This function extracts substring(s) which match a given pattern.
#' @param x Character vector.
#' @param pattern Regular expression.
#' @param sep Character string which separates the fields. 
#' @param mult Logical, if multiple matching fields should be returned,
#' or otherwise NA.
#' @param unlist Logical, unlists multiple results.
#' @param cores Integer for number of computational cores to use.
#' @return
#' A list of character vectors containing the substrings that are
#' matching \code{pattern} and are separated by \code{sep} or \code{NA} if
#' the pattern could not be found.
#' @examples
#' #
#' 
#' s <- c("A1 B1 C1","A2 B2", "AA A1", "AA", "B1 A1", "BB AB A1")
#' strextr(s, "^[AB][[:digit:]]$")
#' strextr(s, "^[AB][[:digit:]]$", mult = TRUE)
#' strextr(s, "^[AB][[:digit:]]$", mult = TRUE, unlist = TRUE)
#' strextr(s, "^[C][[:digit:]]$")
#' 
#' #

#' @export strextr
strextr <- function (x, pattern, sep = " ", mult = F, unlist = F, cores = 1) {  
    x <- strsplit(x, sep)
    x <- mclapply(x, function (y) {
      y <- grep(pattern, y, value = T)
      l <- length(y)
      if ((!mult & l>1) | l==0) y <- NA
      return(y)
    }, mc.cores = cores)
    if (!mult | unlist) x <- unlist(x)
    return(x)
}
