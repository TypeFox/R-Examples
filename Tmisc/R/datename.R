#' Filename with date
#' 
#' Returns string with today's date in YYYY-MM-DD- format concatenated to filename.
#' 
#' @author Stephen Turner
#' @keywords keywords
#' 
#' @param filename A filename string.
#' 
#' @return String with today's date in YYYY-MM-DD- format concatenated to filename.
#' 
#' @examples
#' datename("myfile.png")
#' 
#' @export
datename <- function(filename="filename") paste0(Sys.Date(), "-", filename)
