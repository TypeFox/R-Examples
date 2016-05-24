#' Collect metadata from details call
#'
#' @author Edmund Hart \email{edmund.m.hart@@gmail.com}
#' @param fs_details_obj object
#' @references \url{http://api.figshare.com}
#' @export
#' @import RJSONIO 
#' @examples \dontrun{
#' fs_auth()
#' my_article <- fs_details("138")
#' summary_fs_details(my_article)
#' }
summary_fs_details <- function(fs_details_obj){
  fs_details_vec <- unlist(fs_details_obj)
  fs_summary <- split(fs_details_vec,names(fs_details_vec))
   ### We can change these keys to print what we think is most salient
  print_keys <- c("article_id","defined_type","doi","title","description","shares","views","downloads","owner.full_name","authors.full_name","tags.name","categories.name","files.name","links.link")
  print_names <- c("Article ID", "Article type","DOI","Title","Description","Shares","Views","Downloads","Owner","Authors","Tags","Categories","File names","Links")
  
  for(i in 1:length(print_keys)){
    cat(print_names[i],":",paste(unname(fs_summary[[print_keys[i]]]),sep="",collapse=", "))
    cat("\n")
  }
  
}



