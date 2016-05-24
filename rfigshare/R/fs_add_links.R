#' Add link to article 
#'
#' Adds url links to the metadata of an article 
#' @author Carl Boettiger \email{cboettig@@gmail.com}
#' @param link the url you wish to add (can be list of urls)
#' @param article_id the id number of the article 
#' @param session (optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache as long as authentication has been run. 
#' @param debug logical, should function return details of PUT request?
#' @return output of PUT request (invisibly)
#' @seealso \code{\link{fs_auth}}
#' @references \url{http://api.figshare.com}
#' @import RJSONIO httr
#' @export
#' @examples \dontrun{
#' fs_add_links(138, list("http://carlboettiger.info", "http://ropensci.org")) 
#' }
fs_add_links <- 
function(article_id, link, session = fs_get_auth(), debug = FALSE){
  
  if(is.list(link)){
    link <- unlist(link)
  }
  
  base <- "http://api.figshare.com/v1"
  method <- paste("my_data/articles", article_id, "links", sep= "/")
  request <- paste(base, method, sep = "/")

  for(i in 1:length(link)){
    body <- toJSON(list("link"= link[i]))
    config <- c(config(token = session), 
                add_headers("Content-Type" = "application/json"))
     post <- PUT(request, config = config, body = body)
  }
  if(debug | post$status_code != 200)
    post
  else
    invisible(post)
}



