#' Delete article (private or drafts only) or attached file
#' 
#' @author Carl Boettiger \email{cboettig@@gmail.com}
#' @param article_id the id number of the article 
#' @param file_id the id number of the file, if removing an attached file from a fileset.  
#'   file_id defaults to NULL, removing the entire draft or private article. 
#' @param session (optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache as long as figshare_auth has been run.  
#' @param debug display return value of request?
#' @return output of DELETE request (invisibly) 
#' @seealso \code{\link{fs_auth}}
#' @references \url{http://api.figshare.com}
 
#' @import RJSONIO
#' @export
#' @examples \dontrun{
#' fs_delete(123)
#' 
#' ## Delete all attachments in the second-most-recent entry in my library
#' my_lib <- fs_browse(mine=TRUE)
#' article_id <- my_lib[[2]]$article_id
#' file_ids <- sapply(my_lib[[2]]$files, `[[`, "id")
#' sapply(file_ids, function(id) fs_delete(article_id, id))

#' }
fs_delete <- 
function(article_id, file_id = NULL, session = fs_get_auth(), debug = FALSE){
  base <- "http://api.figshare.com/v1"
  method <- paste("my_data/articles", article_id, sep = "/")
  if(!is.null(file_id))
    method <- paste(method, "files", file_id, sep = "/")
  request <- paste(base, method, sep = "/")
  config <- config(token = session)
  del <- DELETE(request, config = config)
  if(debug | del$status_code != 200)
  invisible(del)
}



