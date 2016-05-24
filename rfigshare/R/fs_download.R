
#' Get details for an article
#'
#' @author Carl Boettiger \email{cboettig@@gmail.com}
#' @param article_id number
#' @param urls_only logical (default TRUE) to only return the URLs to the 
#' downloadable objects but do not call download.file.  If FALSE, will download files
#' @param mine logical (default FALSE). Set to true to see article details for your own non-public articles
#' @param session the authentication credentials from \code{\link{fs_auth}}
#' @param show_versions logical, show what versions are available
#' @param version show a given version number
#' @param ... additional arguments to \code{\link{download.file}}
#' @seealso \code{\link{fs_auth}} \code{\link{download.file}}
#' @references \url{http://api.figshare.com} \url{https://github.com/ropensci/rfigshare}
#' @import httr
#' @export
#' @examples \dontrun{
#' url <- fs_download(90818)
#' data <- read.csv(url)
#' articles <- fs_search("SciFund")
#' ids <- fs_ids(articles)
#' fs_download(ids, urls_only=FALSE)
#' }
fs_download <- 
  function(article_id, urls_only = TRUE, mine=is_mine(article_id), session = fs_get_auth(),
         show_versions=FALSE, version=NULL, ...) {
    details <- lapply(article_id, fs_details, mine = mine, session = session,
           show_versions = show_versions, version = NULL)

    filenames <- unlist(sapply(details, function(output) 
      unlist(lapply(output$files, function(f) f$name))))
    urls <- unlist(sapply(details, function(output) 
      unlist(lapply(output$files, function(f) f$download_url))))

    if(!urls_only)
      sapply(1:length(urls), function(i) 
             download.file(urls[i], destfile=filenames[i], 
                           ..., method = "internal"))
    urls
  }


