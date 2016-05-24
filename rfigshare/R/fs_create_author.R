
#' Creates a figshare author 
#' @author Carl Boettiger \email{cboettig@@gmail.com}
#' @param full_name full name of the author to create
#' @param session (optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache as long as figshare_auth has been run.  
#' @param debug return PUT request visibly?
#' @return author ID numbers
#' @seealso \code{\link{fs_auth}}
#' @references \url{http://api.figshare.com}
#' @import RJSONIO
#' @export
#' @examples \dontrun{
#' fs_create_author("Benjamin Franklin") 
#' } 
fs_create_author <- 
function(full_name, session = fs_get_auth(), debug = FALSE){
  base <- "http://api.figshare.com/v1"
  method <- "my_data/authors"
  request <- paste(base, method, sep="/")
  body <- toJSON(list("full_name"=full_name))
  config <- c(config(token = session), 
            add_headers("Content-Type" = "application/json"))
  request <- build_url(parse_url(request)) # perform % encoding
  post <- POST(request, config = config, body = body)
  if(debug)
    post
  else if(post$status_code == 400){ # Return id if already registered 
    warning(fromJSON(content(post, "text")))
#    fs_author_search(full_name)[[1]]$id # not working at this time
  }
  else if(post$status_code == 201){
    p <- RJSONIO::fromJSON(content(post, "text"))
    p$id
  }
}



