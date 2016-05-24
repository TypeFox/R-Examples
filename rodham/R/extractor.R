#' Get pdf to text extractor (xpdf)
#'
#' @description Downloads and extracts pdf to text extractor, see details.
#'
#' @param dest Destination folder defaults to \code{C:/}
#'
#' @details If function fails you can download the
#' extractor manually from \url{http://www.foolabs.com/xpdf/}
#' (tested on Windows only)
#'
#' @return Returns full path to pdftotext executable
#'
#' @seealso \code{\link{get_emails}}
#'
#' @author John Coene \email{jcoenep@@gmail.com}
#'
#' @export
get_xpdf <- function(dest = "C:/"){
  os <- Sys.info()['sysname'] # get os
  lst <- OStoURI(os) # check os
  temp_zip <- tempfile(fileext = lst$ext) # create temp
  download.file(lst$uri, destfile = temp_zip) # download
  unzip(zipfile = temp_zip, exdir = dest) # unzip
  unlink("temp_zip", recursive=TRUE) # delete temp zip once unzipped
  p <- list.files(paste0(dest, "xpdfbin-win-3.04/bin64/"))
  p <- p[grep("pdftotext", p)]
  message("xpdf successfully uzipped, use: \n",
          dest, "xpdfbin-win-3.04/bin64/", p, "\n",
          "as extractor in get_emails")
  return(paste0(dest, "xpdfbin-win-3.04/bin64/", p))
}

