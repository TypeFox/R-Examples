#' Make an article public (for private or draft articles)
#'
#' @author Carl Boettiger \email{cboettig@@gmail.com}
#' @param article_id the id number of the article 
#' @param session (optional) the authentication credentials from \code{\link{fs_auth}}. 
#' @return output of PUT request (invisibly)
#' @details NOTE: Public articles are assigned DOIs and cannot be deleted or made private once declared public! Public articles do not count against your quota space.  
#' @seealso \code{\link{fs_auth}}
#' @references \url{http://api.figshare.com}
#' @export
#' @examples \dontrun{
#' fs_make_public(123)
#' }
fs_make_public <- function(article_id, session = fs_get_auth()){
  base <- "http://api.figshare.com/v1"
  method <- paste("my_data/articles", article_id, "action/make_public", sep="/")
    request = paste(base, method, sep="/")
  POST(request, config(token = session))
}


