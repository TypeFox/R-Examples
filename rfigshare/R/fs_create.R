#' Create a FigShare article (draft)
#' 
#' Articles must be created with \code{\link{fs_create}}
#' with essential metadata.  Then you can add files with
#' \code{\link{fs_upload}}, add categories, tags or authors
#' with \code{\link{fs_add_categories}} or \code{\link{fs_add_tags}}
#' \code{\link{fs_add_authors}}.  Authors not registered with a FigShare
#' id can be created with \code{\link{fs_create_author}}.  You can
#' edit the original metadata with \code{\link{fs_update}}.
#' Finally, release the article as either private or public with
#' \code{\link{fs_make_private}} or \code{\link{fs_make_public}}.
#' Before creating the article, you must authenticate using
#' \code{\link{fs_auth}}.
#' @param title for the article
#' @param description of the article
#' @param type one of: dataset, figure, media, poster, paper or fileset. (Only filesets can have multiple uploaded files attached).  
#' @param session the authentication credentials from \code{\link{fs_auth}}
#' @param debug print full post call return
#' @return article id 
#' @seealso \code{\link{fs_auth}}, \code{\link{fs_upload}}
#' @references \url{http://api.figshare.com}
#' @import RJSONIO
#' @import methods
#' @import httr
#' @export
#' @examples \dontrun{
#' fs_create("My Title", "A description of the object", "dataset")
#' }
fs_create <- 
function(title, description, type = 
         c("dataset", "figure", "media", "poster", "paper", "fileset"),
         session = fs_get_auth(), debug = FALSE) {
# TODO: Return (or at least message) the article ID number.  Error handling for types?
  type <- match.arg(type)
  base <- "http://api.figshare.com/v1"
  method <- "my_data/articles"
  request <- paste(base, method, sep="/")
  meta <- toJSON(list("title" = title, 
                      "description" = description, 
                      "defined_type" = type))
  config <- c(config(token = session), 
              add_headers("Content-Type" = "application/json"))
  post <- POST(request, config = config, body = meta)
  if(debug | post$status_code != 200)
    post
  else {
    p <- RJSONIO::fromJSON(content(post, "text"))
    article_id <- p$article_id

    if(is.numeric(article_id))
      message(paste("Your article has been created! Your id number is", article_id))

    article_id
  }
}



