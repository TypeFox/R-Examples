#' Attempts to download a PDF using a DOI link.
#'
#' Tries to download a PDF file using the digital objected identifier (DOI) link.
#' Uses ad hoc searches of journal HTML pages to detect candidate PDFs for 
#' download.  Downloads all candidate pdfs.
#'
#' @param DOI A string of the DOI (digital object identifier) used to identify
#'    the source of a journal article PDF file(s).
#' @param directory A string of the location (directory) were downloaded PDF
#'    files are to be saved.  Directory name must end with "\\\\". 
#' @param theFileName Used to rename the downloaded file.  No need to include
#'    extension ".pdf".
#' @param validatePDF When \code{"TRUE"} will only save to files that are valid
#'    PDF documents.  When \code{"FALSE"} will save all candidate files, even if
#'    they are not valid PDF formats.
#' @param quiet When \code{"FALSE"} does not print to console download progress 
#'    and summary.
#'
#' @return A string describing the download success.  If unsuccessful,
#'    returns the type of error during the download attempt.
#'
#' @seealso \code{\link{PDFs_collect}}
#'
#' @importFrom stringr str_extract
#' @importFrom utils download.file
#' @export PDF_download

PDF_download <- function(DOI, 
                         directory = getwd(), 
                         theFileName = "temp", 
                         validatePDF = TRUE, 
                         quiet = FALSE) {
  
  
  if(!quiet) {
    message(paste0("Collecting PDF from DOI: ", DOI))
    message(paste0("\t\t\tExtraction 1 of 2: HTML script...."), appendLF = FALSE) 
  }
  
  if(is.URLconnectable(paste0("http://dx.doi.org/", DOI))) {
    urlMessage <- " successful"
    theHTMLvector <- getHTMLfromURL(paste0("http://dx.doi.org/", DOI))
    
    if(!quiet) {
      message(paste0(urlMessage)) 
      message(paste0("\t\t\tExtraction 2 of 2: PDF download..."), appendLF = FALSE) 
    }
    
    wasPDFdownloaded <- extractPDFsFromHTML(theHTMLvector, 
                                            directory, 
                                            theFileName, 
                                            validatePDF)
    
    if(wasPDFdownloaded == TRUE) {
      downloadMessage <- " successful"
      downloadOutcome <- "downloaded"
    } else {
      downloadMessage <- wasPDFdownloaded
      downloadOutcome <- "download error"
    }
    
     if(!quiet) message(paste0(downloadMessage, 
                               ifelse(downloadOutcome == " downloaded", 
                                      paste0(" (filename: ", theFileName, ".pdf)"), "")))
    
  
  } else {
    urlMessage <- " cannot open: HTTP status was '404 Not Found'"
    downloadMessage <- " skipped"
    
    if(!quiet) {
      message(paste0(urlMessage)) 
      message(paste0("\t\t\tExtraction 2 of 2: PDF download...", 
                     downloadMessage))
    }
    
    downloadOutcome <- "URL error"
  }
    
  return(downloadOutcome)
}