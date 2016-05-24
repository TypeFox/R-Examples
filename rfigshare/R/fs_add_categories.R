#' Add a category to article
#' 
#' @author Edmund Hart \email{edmund.m.hart@@gmail.com}
#' @param article_id the id number of the article 
#' @param category_id is a vector of integers corresponding to categories or a vector of category names  
#' @param session (optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache as long as figshare_auth has been run.  
#' @param debug return PUT results visibly?
#' @return output of PUT request (invisibly)
#' @seealso \code{\link{fs_auth}}
#' @references \url{http://api.figshare.com}
#' @import RJSONIO httr
#' @export
#' @examples \dontrun{
#' fs_add_categories(138, "Ecology")
#' }
fs_add_categories <- function(article_id, category_id, 
                              session = fs_get_auth(), 
                              debug = FALSE){
  
  if(is.list(category_id)){
    category_id <- unlist(category_id)
  }
  suppressWarnings(my_ids <- as.numeric(category_id))
  if(sum(is.na(my_ids)) > 0){
    category_id <- fs_cat_to_id(category_id)
  }
  
  base <- "http://api.figshare.com/v1"
  method <- paste("my_data/articles", article_id, "categories", sep = "/")
  request <- paste(base, method, sep = "/")
  
  for(i in 1:length(category_id)){
    body <- toJSON(list("category_id" = category_id[i]))
    config <- c(config(token = session), 
                add_headers("Content-Type" = "application/json"))
    
    post <- PUT(request, config = config, body = body)
    if(debug | post$status_code != 200)
      post
    else
      invisible(post)
  }
}



#'  Helper function that matches string categories to id's
#' 
#' @author Edmund Hart \email{edmund.m.hart@@gmail.com}
#' @param category_id Must be a valid category string, regardless of case 
#' @return a vector of integers corresponding to valid figshare categories
#' @references \url{http://api.figshare.com}
#' @import RJSONIO httr plyr
fs_cat_to_id <- function(category_id){
  if(!exists("cat_names")) 
    cat_names <- RJSONIO::fromJSON(content(GET("http://api.figshare.com/v1/categories"), "text"))
  name_db <- ldply(cat_names$items, data.frame)  
  name_db$name <- tolower(name_db$name)
  my_matches <- match(tolower(category_id), name_db$name)
  if(sum(is.na(my_matches)) > 0){
    stop("One or more of your input strings is not a valid category") 
  }
  if(sum(is.na(my_matches)) == 0){
    return(name_db$id[my_matches]) 
  }
}


