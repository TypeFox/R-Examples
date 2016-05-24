effort_save <- function (aDataFrame, 
                         column_name = "REVIEWERS", 
                         directory = getwd(), 
                         quiet = FALSE) {
  
  split_List <- split(aDataFrame, f = aDataFrame[, column_name])
  lapply(split_List, function(x, column_name) {
    fileName <- paste0("effort_", x[1, column_name], ".csv")
    if(file.exists(file.path(directory, fileName))) {
      while (file.exists(file.path(directory, fileName)) != FALSE) fileName <- renameFile(fileName)
      if(!quiet) message(paste0("File already existed, renamed: ", fileName))
    }
    fileName <- file.path(directory, fileName)
    write.csv(x, file = fileName, na = "NA", row.names = FALSE)
  }, column_name = column_name)
  
  if(!quiet) message(paste0(length(split_List), " files saved in: ", directory))
}

Landis_Koch_scale <- function(aKappa) {
  if(aKappa <= 0) {
    return("poor")
  } else if (aKappa <= 0.2) {
    return("slight")
  } else if (aKappa <= 0.4) {
    return("fair")
  } else if (aKappa <= 0.6) {
    return("moderate")
  } else if (aKappa <= 0.8) {
    return("substantial")
  } else return("almost perfect")
}