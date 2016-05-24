#' Add a tag to an article
#' 
#' @author Carl Boettiger \email{cboettig@@gmail.com}
#' @param article_id the id number of the article to create
#' @param tag name of the tag to add (or list of tags)
#' @param session the authentication credentials from \code{\link{fs_auth}}
#' @return output of PUT request (invisibly)
#' @seealso \code{\link{fs_auth}}
#' @references \url{http://api.figshare.com}
#' @import RJSONIO
#' @export
#' @keywords internal
#' @examples \dontrun{
#'  fs_add_tag(138, "phylogenetics") 
#' }
fs_add_tags <- 
  function(article_id, tag, session = fs_get_auth(), debug = FALSE){
    if(is.list(tag)){
      tag <- unlist(tag)
    }
    base <- "http://api.figshare.com/v1"
    method <- paste("my_data/articles", article_id, "tags", sep = "/")
    request <- paste(base, method, sep="/")
    for(i in 1:length(tag)){
      body <- toJSON(list("tag_name" = tag[i]))
      config <- c(config(token = session), 
                  add_headers("Content-Type" = "application/json"))
      
      post <- PUT(request, config = config, body = body)
      if(debug | post$status_code != 200)
        post
      else
      invisible(post)
    }
  #  fs_tag_as_rfigshare(article_id)
  }
