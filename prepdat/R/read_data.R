#' Reads a File in a txt or csv Format that Contains a Table and Creates a Data
#' Frame from it
#' 
#' @param file_name A string with the name of the file to be read into R.
#' @param notification Logical. If TRUE, prints messages about the progress of
#'   the function. Default is \code{TRUE}.
#' @return A data frame of the table specified in \code{file_name}.
read_data <- function(file_name, notification = TRUE) {
  
  # Get file_name extension 
  extension <- substr(file_name, nchar(file_name) - 3, nchar(file_name))
  
  # Check if extension is ".txt" or ".csv"
  if (extension == ".txt") {
    if (file_name %in% list.files()) {
      if (notification == TRUE) {
        # Message
        message(paste("Reading", file_name, "into R"))
      }
      # Return
      return(read.table(file_name, header = TRUE))
    } else {
      # Error handling 
      stop(paste("Oops! Did not find", file_name, "in", getwd()))
    }
  } else if (extension == ".csv") {
    if (file_name %in% list.files()) {
      if (notification == TRUE) {
        # Message
        message(paste("Reading", file_name, "file into R"))
      }
      # Return
      return(read.csv(file_name, header = TRUE))
    } else {
      # Error handling 
      stop(paste("Oops! Did not find", file_name, "in", getwd()))
    }
  } else {
    # Error handling
    stop("Oops! file_name should end with txt or csv extension (and be in that format)")     
  }
}