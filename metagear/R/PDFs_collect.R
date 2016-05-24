#' Attempts to download PDFs from multiple DOI links.
#'
#' Tries to download a collection of PDF files using multiple digital object 
#' identifier (DOI) links.  Updates a data frame with the success of these 
#' downloads.  The function is a wrapper for \code{\link{PDF_download}}.  NOTE: A 
#' single DOI may generate multiple PDF files.
#'
#' @param aDataFrame A data frame containing a column of DOIs and a column of
#'    individual file names for each downloaded PDF.
#' @param DOIcolumn The label of the column containing all the DOI links.
#' @param FileNamecolumn The label of the column containing all the strings 
#'    that will be used to rename the downloaded files.
#' @param directory A string of the location (directory) were downloaded PDF
#'    files are to be saved.  NOTE: helps to have this directory created before
#'    initializing the \code{PDFs_collect} function.
#' @param validatePDF When \code{TRUE} will only save to files that are valid
#'    PDF documents.  When \code{FALSE} will save all candidate files, even if
#'    they are not valid PDF formats.
#' @param quiet When \code{FALSE} does not print to console individual 
#'    download progress and summary.
#' @param showSummary When \code{FALSE} does not print overall summary of download
#'    successes and failures.
#'
#' @return The data frame with new column containing download-outcome successes.
#'
#' @examples \dontrun{
#'
#' data(example_references_metagear)
#' someRefs <- effort_initialize(example_references_metagear)  
#' dir.create("metagear_downloads")      
#' PDFs_collect(aDataFrame = someRefs, DOIcolumn = "DOI", 
#'              FileNamecolumn = "STUDY_ID", directory = "metagear_downloads")
#' }
#'
#' @seealso \code{\link{PDF_download}}
#'
#' @importFrom utils download.file 
#' @export PDFs_collect

PDFs_collect <- function(aDataFrame, 
                         DOIcolumn, 
                         FileNamecolumn, 
                         directory = "downloads", 
                         validatePDF = TRUE, 
                         quiet = FALSE, 
                         showSummary = TRUE) {

  downloadOutcomes <- mapply(
    function(DOI, directory, theFileName) PDF_download(DOI, 
                                                       directory, 
                                                       theFileName, 
                                                       validatePDF, 
                                                       quiet), 
    aDataFrame[, DOIcolumn], 
    directory, 
    aDataFrame[, FileNamecolumn]
  )
  
  if(showSummary & !quiet) {
    message("\nPDF download summary")
    a <- summary(as.factor(downloadOutcomes))
    for(i in 1:length(a)) message(paste("\t", a[i], "=", names(a[i])))
    #need to fix, to get correct path
    #message(paste0("\tDownloads located in: ", as.character(getwd()), "/", directory))
    message(paste0("\tDownloads located in: ", directory))
  }
  
  aDataFrame$downloadOutcomes <- downloadOutcomes
  return(aDataFrame)
}

# NOTE: tried to parallelize downloads but multiple simultaneous 
# connections failed (or perhaps blocked by  firewall) <- 11/10/2014