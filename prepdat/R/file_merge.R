#' Vertically Merge Files in a Directory into a Single Large Dataset
#'
#' @description Vertically concatenates files containing data tables in a long
#'   format into a single large dataset. In order for the function to work, all
#'   files you wish to merge should be put in one folder containing nothing but
#'   these files. In addition, the working directory should be set to that
#'   folder and all files should be in the same format (either txt or csv).
#'   This function is very useful for concatenating raw data files of
#'   individual subjects in an experiment (in which each line corresponds to a
#'   single observation in the experiment) to one raw data file that includes
#'   all subjects.
#' @param file_name A string with the name of the file of the large dataset the
#'   function returns in case \code{save_table} is \code{TRUE}. Default is
#'   \code{"dataset.txt"}.
#' @param save_table Logical. If TRUE, saves the dataset into a txt file.
#'   Default is \code{"dataset.txt"}.
#' @param notification Logical. If TRUE, prints messages about the progress of
#'   the function. Default is \code{TRUE}.
#' @return The merged dataset \code{dataset}.
file_merge <- function(file_name = "dataset.txt", save_table = TRUE,
                       notification = TRUE) {

  # Error handling
  # Get file_name extension
  extension <- substr(file_name, nchar(file_name) - 3, nchar(file_name))
  if (extension != ".txt") {
    stop("Oops! file_name must end with txt extension")
  }

  # Make a list of all files in the directory
  file_list <- list.files()

  # Message how many files were found
  if (notification == TRUE) {
    if (length(file_list) == 1) {
      message(paste("Found", length(file_list), "file"))
    } else {
      message(paste("Found", length(file_list), "files"))
    }
  }

  # Counter for loop
  i <- 1

  # Merge files vertically
  for (file in file_list) {

    # Get file extension
    extension <- substr(file, nchar(file) - 3, nchar(file))

    # Read first file in file_list into dataset
    if (match(file_list[i], file_list) == 1) {
      if (extension == ".txt") {
        dataset <- read.table(file, header = TRUE)
      } else if (extension == ".csv") {
        dataset <- read.csv(file, header = TRUE)
      }
    }

    # Append current file to large dataset if it is not the first file
    if (match(file_list[i], file_list) != 1) {
      # Read current file into temp_dataset
      if (extension == ".txt") {
        temp_dataset <- read.table(file, header = TRUE)
      } else if (extension == ".csv") {
        temp_dataset <- read.csv(file, header = TRUE)
      }
      # Append temp_dataset to dataset
      dataset <- rbind(dataset, temp_dataset)
      # Remove temp_dataset
      rm(temp_dataset)
    }

    # Increase counter
    i <- i + 1
  }  # End of for loop

  # Save table in case save_table is set to TRUE
  if (save_table == TRUE) {
    if (notification == TRUE) {
      # Message
      message(paste(length(list.files()), "were merged and saved into", file_name))
    }
    write.table(dataset, row.names = FALSE, file = file_name)
  }

  if (notification == TRUE) {
    # Message
    message("file_merge() finished!")
  }

  # Return
  return(dataset)
}
