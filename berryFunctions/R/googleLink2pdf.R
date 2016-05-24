#' extract pdf link from google search result
#' 
#' restrict pdf link from a google search to actual link with text processing
#' 
#' @return Characterstring with only the basic link
#' @note The function is not vectorized! If you have many links, use a loop around this function... 
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2012
#' @seealso \code{\link{strsplit}}, \code{\link{gsub}}
#' @keywords character
#' @export
#' @examples
#' 
#' Link <- paste0("http://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1",
#'         "&cad=rja&sqi=2&ved=0CDIQFjAA&url=http%3A%2F%2Fcran.r-project.org",
#'         "%2Fdoc%2Fmanuals%2FR-intro.pdf&ei=Nyl4UfHeOIXCswa6pIC4CA",
#'         "&usg=AFQjCNGejDwPlor4togQZmQEQv72cK9z8A&bvm=bv.45580626,d.Yms")
#' googleLink2pdf(Link)
#' 
#' @param googlelink Character string: A search result address
#' 
googleLink2pdf <- function(
googlelink
)
{
pdflink <- strsplit(googlelink, "&url=")[[1]][2]
pdflink <- strsplit(pdflink, "&ei=")[[1]][1]
pdflink <- gsub("%24", "$", pdflink)
pdflink <- gsub("%2F", "/", pdflink)
pdflink <- gsub("%3A", ":", pdflink)
pdflink <- gsub("%25", "%", pdflink)
pdflink
}
