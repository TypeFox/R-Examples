#' wget download a file from internet and save it locally
#' @title Download a file from internet and save it locally
#' @author Marc Girondot
#' @return Nothing
#' @param url The url where to download file
#' @param ... The parameters send to download.file()
#' @description Download a file from internet and save it locally. This function is a wrapper for
#' download.files() that keep the name identical and can get several files at once.
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' # Save locally the files send in the parameter url
#' wget(c("https://cran.r-project.org/web/packages/HelpersMG/HelpersMG.pdf", 
#'          "https://cran.r-project.org/web/packages/embryogrowth/embryogrowth.pdf"))
#' }
#' @export


wget <- function(url=stop("At least one internet adress is required"), ...) {
  for (i in 1:length(url))
    do.call(download.file, modifyList(list(url=url[i], destfile=basename(url[i])), list(...)))
  }
