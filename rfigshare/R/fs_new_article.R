

#' Create a FigShare article.  
#' 
#' fs_new_article is a wrapper around many other rfigshare commands to provide convenient posting.   
#' @param title for the article, see \code{\link{fs_create}} for details.
#' @param description of the article, see \code{\link{fs_create}} for details.
#' @param type one of: dataset, figure, media, poster, or paper, see \code{\link{fs_create}} for details. 
#' @param authors Orded list of authors for the article, see \code{\link{fs_add_authors}} for details
#' @param categories list of categories or category id numbers, see \code{\link{fs_add_categories}} for details.   
#' @param tags list of tags, see \code{\link{fs_add_tags}} for details.   
#' @param links list of links to add, see \code{\link{fs_add_links}} for details
#' @param files path to the files to add, see \code{\link{fs_upload}} for details
#' @param visibility one of "draft", "private" or "public".  A draft document can still be edited and modified.  A public document is visible to everyone and cannot be deleted (though additional authors to the work can still "claim" their authorship).  
#' @param session (optional) credentials, see \code{link{fs_auth}}
#' @return article id 
#' @seealso \code{\link{fs_auth}}, \code{\link{fs_add_categories}}, \code{\link{fs_add_authors}}, \code{\link{fs_add_tags}}, \code{\link{fs_add_links}}
#' @references \url{http://api.figshare.com}
#' @import RJSONIO
#' @export
#' @examples \dontrun{
#' write.csv(mtcars, "mtcars.csv")
#' id <- fs_new_article(title="A Test of rfigshare", 
#'                     description="This is a test of the fs_new_article function and related methods", 
#'                     type="dataset", 
#'                     authors=c("Karthik Ram", "Scott Chamberlain"), 
#'                     tags=c("ecology", "openscience"), 
#'                     categories="Ecology", 
#'                     links="http://ropensci.org", 
#'                     files="mtcars.csv",
#'                     visibility="private")
#' }

fs_new_article <- function(title, description, type = 
         c("dataset", "figure", "media", "poster", "paper", "fileset"),
         authors = NA, categories = NA, tags = NA, 
         links = NA, files = NA, 
         visibility = c("draft", "private", "public"),
         session = fs_get_auth()){

  article_id <- fs_create(title = title, 
                          description = description, 
                          type = type, 
                          session = session)
  visibility <- match.arg(visibility)
  if(!is.nas(authors))
    fs_add_authors(article_id, authors, session)
  if(!is.nas(categories))
    fs_add_categories(article_id, categories, session)
  if(!is.nas(tags))
    fs_add_tags(article_id, tags, session)
  if(!is.nas(links))
    fs_add_links(article_id, links, session)
  if(!is.nas(files))
    fs_upload(article_id, files, session)
  if(visibility == "private")
    fs_make_private(article_id, session)
  if(visibility == "public")
    fs_make_public(article_id, session)

#   fs_tag_as_rfigshare(article_id) # performed by fs_upload
  article_id
}

is.nas <- function(x) any(sapply(x, is.na))

