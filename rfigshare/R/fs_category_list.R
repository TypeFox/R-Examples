
#' List all categories
#'
#' @return a table of all the categories 
#' @references \url{http://api.figshare.com}
#' @import httr
#' @param debug enable debugging
#' @export
#' @examples \dontrun{
#' fs_categories_list()
#' }
fs_category_list <- function(debug = FALSE){
    response <- GET("http://api.figshare.com/v1/categories")
    if(debug | response$status_code != 200)
      response
    else{ 
      a <- fromJSON(content(response, "text"))
      b <- t(sapply(a$items, function(d) c(id = d$id, name = d$name)))
      b <- as.table(b)
      rownames(b) <- NULL
      b
    }
}
