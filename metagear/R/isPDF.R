#' Evaluates whether a file is a PDF document.
#'
#' Checks if provided file is in Portable Document Format (PDF).    
#'
#' @param aFileName A string that identifies a file name (and directory path) of
#' the PDF candidate.  
#'
#' @return A logical value indicating whether the file is a PDF document.
#'
#' @importFrom stringr str_extract
#' @export isPDF

isPDF <- function(aFileName) {
  fileContents <- readLines(aFileName, n = -1, ok = TRUE, warn = FALSE)
  return(any(!is.na(str_extract(fileContents, ".*%PDF-1.*"))))
}