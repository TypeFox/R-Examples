#' Find WWW Domains
#'
#' @description   Hand function to return the \code{www} domain.
#'
#' @param x A vector from which the domain information is desired.
#'
#' @keywords Manipulation
#'
#' @export
#' @examples
#'  x1 <- "http://stackoverflow.com/questions/19020749/function-to-extract-domain-name-from-url-in-r"
#' x2 <- "http://www.talkstats.com/"
#' x3 <- "www.google.com"
#'
#' has.domain(x3)
#'
#' sapply(list(x1, x2, x3), has.domain)

`has.domain` <- function(x){
  x<-tolower(x)
 out <- strsplit(gsub("http://|https://|www\\.", "", x), "/")[[c(1, 1)]]
return(out)
 }
