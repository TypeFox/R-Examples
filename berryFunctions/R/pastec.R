#' Paste with collapse = ", "
#' 
#' Helper function \code{\link{paste}} with collapse = ", "
#' 
#' @return Single character string
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, April 2015
#' @seealso \code{\link{paste}}, \code{\link{paste0}}
#' @keywords character
#' @export
#' @examples
#' 
#' listoferrors <- c("filetype", "header", "nonemptyline")
#' message("The following entities were corrupted:\n", pastec(listoferrors))
#' pastec("Part1", c("Part2", "Part3"), letters[1:3])
#' 
#' @param \dots Object(s) to be \code{\link{paste}d} to a character vector
#' @param sep Character string to separate single strings. DEFAULT: " "
#' @param collapse Character string between combined strings. DEFAULT: ", "
#' 
pastec <- function(
...,
sep=" ",
collapse=", "
)
{
paste(c(...), sep=sep, collapse=collapse)
}

## Old idea to enable newline insertion every n elements, apparently not necessary.
##newlines=TRUE # try to add line breaks at suitable intervals
##\item{newlines}{Add line breaks? (ca. at console width). DEFAULT: TRUE}
##dd <- paste(..., sep=sep, collapse=collapse)
##if(newlines)
##{
##  w <- getOption("width")
##  if(nchar(dd) < w) return(dd) 
##  dd2 <- strsplit(dd, collapse)[[1]]
##  sel <- cumsum(nchar(dd2)+2) < w
##  dd3 <- dd2[sel]
##}
##return(dd)
