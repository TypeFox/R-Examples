#' Get a list of article id numbers from a search return 
#' 
#' @author Carl Boettiger \email{cboettig@@gmail.com}
#' @param object the output of a search 
#' @return a list of article id numbers 
#' @references \url{http://api.figshare.com}
#' @import httr
#' @export
#' @examples \dontrun{
#' figshare_category() 
#' }
fs_ids <- function(object){
#  a <- parsed_content(object)
#  sapply(a$items, `[[`, "article_id")
  sapply(object, `[[`,  "article_id")
}
