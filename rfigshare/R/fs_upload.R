#'  Upload file to an article 
#' 
#' @details Article must be a draft, i.e. created by \code{\link{fs_create}} and not yet made public or private. Only articles of type "fileset" can have multiple files uploaded.  
#' @author Carl Boettiger \email{cboettig@@gmail.com}
#' @param article_id an article id number or a character string (or list) of numbers
#' @param file path to file to upload, or character string (or list) of files (paths)
#' @param session the authentication credentials from \code{\link{fs_auth}} (optional)
#' @details If only a single id number is given but a character string of files is given,
#' then be sure that the id corresponds to an object of type "fileset".  If article_id list
#' has more than one id, then there must be a corresponding file path for each id.  
#' @seealso \code{\link{fs_auth}}
#' @references \url{http://api.figshare.com}
#' @examples
#' \dontrun{
#' id <- fs_create("Title", "description", "figure")
#' fs_upload(id, "file.png")
#' } 
#' @export
fs_upload <- function(article_id, file, session = fs_get_auth()) {
  # handle lists by converting to expected type.  
  if(is(article_id, "list"))
    article_id <- as.character(article_id)
  if(is(file, "list"))
    file <- as.character(file)

  if(length(article_id) == 1 && length(file) == 1)
    fs_upload_one(article_id, file, session)
  if(length(article_id) > 1)
    sapply(1:length(article_id), function(i) fs_upload_one(article_id[i], file[i], session))
  if(length(article_id) == 1 && length(file) > 1)
    sapply(file, function(f) fs_upload_one(article_id, f, session))


}




#'  Upload file to an article 
#' 
#' @details Article must be a draft, i.e. created by \code{\link{fs_create}} and not yet made public or private. Only articles of type "fileset" can have multiple files uploaded.  
#' @author Carl Boettiger \email{cboettig@@gmail.com}
#' @param article_id number
#' @param file path to file to upload
#' @param session the authentication credentials from \code{\link{fs_auth}} (optional)
#' @seealso \code{\link{fs_auth}}
#' @references \url{http://api.figshare.com}
#' @keywords internal
#' @examples
#' \dontrun{
#' id <- fs_create("Title", "description", "figure")
#' fs_upload(id, "file.png")
#' } 
fs_upload_one  <- function(article_id, file, session = fs_get_auth()) {
  base <- "http://api.figshare.com/v1"
  method <- paste("/my_data/articles", article_id, "files", sep = "/")
  request <- paste(base, method, sep="")
  body <- list(filedata = upload_file(file))
  config <- c(config(token = session), add_headers("Content-Type" = "multipart/form-data"))
  out <- PUT(request, config = config, body = body)
  fs_tag_as_rfigshare(article_id)
  out
}
