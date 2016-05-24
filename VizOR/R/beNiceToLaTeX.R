##' Sanitize deparsed R expressions for inclusion in LaTeX documents.
##' Percent signs and ampersands appearing in R expressions must be escaped for
##' inclusion in LaTeX documents.
##' 
##' @param str Typically, a deparsed R expression to be included in, e.g., a
##'   figure caption in a LaTeX document.
##' @return A string is returned, with LaTeX-unfriendly characters
##'   backslash-escaped.
##' @note This function is a candidate for documentation as package-internal.
##' @author David C. Norris
##' @keywords internal
##'
beNiceToLaTeX <-
function(str){ # Escape LaTeX-unfriendly characters
  str <- gsub("[&]", "\\\\&", str)
  str <- gsub("[%]", "\\\\%", str)
  str
}

