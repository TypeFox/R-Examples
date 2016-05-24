#' Update article title, description, or type 
#'
#' Updates the article title, description or type. If any is not specified, it will remain unchanged.  
#' @author Carl Boettiger \email{cboettig@@gmail.com}
#' @param article_id the id number of the article 
#' @param title for the article (to replace original title)
#' @param mine Set to \code{TRUE} if it refers to an item on your own account
#' @param description of the article (replaces original designation)
#' @param type one of: dataset, figure, media, poster, or paper (replaces original designation)
#' @param session (optional) the authentication credentials from \code{\link{fs_auth}}. 
#' @param debug return httr PUT request visibly?
#' @return output of PUT request (invisibly)
#' @seealso \code{\link{fs_auth}}, \code{\link{fs_add_tags}}
#' @references \url{http://api.figshare.com}
#' @details Updates the title, description, and type of an article.  
#' @import RJSONIO
#' @export
#' @examples \dontrun{
#' fs_update(138, title = "New title") 
#' }
fs_update <- 
function(article_id, title = NA, description = NA, type  =  NA, mine = TRUE,
         session = fs_get_auth(), debug = FALSE){

  ## grab the article details and use those as defaults
  details <- fs_details(article_id, mine = mine, 
                        session = session)
  if(is.na(title))
    title <- details$title
  if(is.na(description))
    description <- details$description
  if(is.na(type))
    type <- details$defined_type

  base <- "http://api.figshare.com/v1"
  method <- paste("my_data/articles", article_id, sep = "/")
  request <- paste(base, method, sep = "/")
  meta <- toJSON(list("title" = title, 
                      "description" = description, 
                      "defined_type" = type))
  config <- c(config(token = session), 
              add_headers("Content-Type" = "application/json"))
  post <- PUT(request, config = config, body = meta)
  if(debug)
    post
  else
    invisible(post)
}



