################################################################################
# TODO LIST
# TODO: Update to use 'fread' when problem with colClasses is solved.


################################################################################
# CHANGE LOG (last 20 changes)
# 09.01.2016: Added more attributes to result.
# 15.12.2015: Removed "0" from the default 'na.strings'.
# 04.12.2015: Added parameter 'na.strings'.
# 09.11.2015: Added "0" to 'na.strings' in 'read.table'.
# 06.10.2015: Added call to 'colConvert' to convert known numeric columns.
# 05.10.2015: Added attributes.
# 31.08.2015: Removed option to manually pick folder using 'choose.dir'.
# 29.08.2015: Added importFrom.
# 01.06.2015: Re-named column name 'File' to 'File.Name' to increase specificity in 'trim'.
# 23.05.2015: Changed names on parameters 'file.name' -> 'import.file'.
# 23.05.2015: Added parameters for auto trim and auto slim.
# 22.05.2015: Added parameters 'file.name', 'time.stamp', and 'ignore.case'.
# 22.05.2015: Re-wrote import loop.
# 15.12.2014: Changed parameter names to format: lower.case
# 20.01.2014: Added parameter 'colClasses = "character"' to 'read.table'.
# 15.01.2014: Added message to show progress.
# 13.01.2014: Added parameter 'na.strings = c("NA","")' to 'read.table'.
# 13.01.2014: Fixed bug when no matching files in folder.
# 10.12.2013: Changed names on parameters 'resultFiles' -> 'file.name'
#              and 'resultFolder' -> 'folder.name'.
# 12.11.2013: Changed 'rbind' to 'rbind.fill' from package 'plyr'.
# 13.06.2013: Added parameter 'debug'. Fixed regexbug when importing from folder.

#' @title Import Data
#'
#' @description
#' Import text files and apply post processing.
#'
#' @details
#' Imports text files (e.g. GeneMapper results exported as text files)
#' as data frames. Options to import one or multiple files. For multiple
#' files it is possible to specify prefix, suffix, and file extension
#' to create a file name filter. The file name and/or file time stamp
#' can be imported.
#' NB! Empty strings ("") and NA strings ("NA") are converted to NA.
#' See \code{\link{list.files}} and \code{\link{read.table}} for additional details.
#' 
#' @param folder logical, TRUE all files in folder will be imported,
#' FALSE only selected file will be imported.
#' @param suffix string, only files with specified suffix will be imported.
#' @param prefix string, only files with specified prefix will be imported.
#' @param import.file string if file name is provided file will be imported
#' without showing the file open dialogue. 
#' @param folder.name string if folder name is provided files in folder
#' will be imported without showing the select folder dialogue. 
#' @param extension string providing the file extension.
#' @param file.name logical if TRUE the file name is written in a column 'File.Name'.
#' NB! Any existing 'File.Name' column is overwritten.
#' @param time.stamp logical if TRUE the file modified time stamp is written
#' in a column 'Time'.
#' NB! Any existing 'Time' column is overwritten.
#' @param separator character for the delimiter used to separate columns
#'  (see 'sep' in \code{\link{read.table}} for details).
#' @param ignore.case logical indicating if case should be ignored. Only applies
#' to multiple file import option.
#' @param auto.trim logical indicating if dataset should be trimmed.
#' @param trim.samples character vector with sample names to trim.
#' @param trim.invert logical indicating if samples should be keept (TRUE) or
#'  removed (FALSE).
#' @param auto.slim logical indicating if dataset should be slimmed.
#' @param slim.na logical indicating if rows without data should remain.
#' @param na.strings character vector with strings to be replaced by NA.
#' @param debug logical indicating printing debug information.
#' 
#' 
#' @return data.frame with imported result.
#' 
#' @export
#' 
#' @importFrom plyr rbind.fill
#' @importFrom utils read.table
# @importFrom data.table fread
#' 
#' @seealso \code{\link{trim}}, \code{\link{slim}}, \code{\link{list.files}}, \code{\link{read.table}}



import <- function (folder = TRUE, extension="txt", 
                    suffix = NA, prefix = NA, 
                    import.file = NA, folder.name = NA,
                    file.name = TRUE, time.stamp = TRUE,
                    separator = "\t", ignore.case = TRUE,
                    auto.trim = FALSE, trim.samples = NULL,
                    trim.invert = FALSE,
                    auto.slim = FALSE, slim.na = TRUE,
                    na.strings = c("NA",""),
                    debug = FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  if(debug){
    print("folder")
    print(folder)
    print("extension")
    print(extension)
    print("suffix")
    print(suffix)
    print("prefix")
    print(prefix)
    print("import.file")
    print(import.file)
    print("folder.name")
    print(folder.name)
    print("file.name")
    print(file.name)
    print("time.stamp")
    print(time.stamp)
    print("ignore.case")
    print(ignore.case)
  }
  
  # Check data ----------------------------------------------------------------
  
  # Check data.
  if(is.na(import.file) && is.na(folder.name)){
    stop("Either 'import.file' or 'folder.name' must be provided")
  }
  
  # Import --------------------------------------------------------------------
  
  # Initialise result data.frame (no match return empty dataframe)
  res <- data.frame()
  
  # Check if result files in folder.
  if (folder) {

    # Get path.    
    folder <- folder.name

    # Check if folder is specified.
    if (!is.na(folder)) {
      
      # Create file filter.
      fileFilter <- paste(".*", sep="")
      if (!is.na(prefix) && nchar(prefix) > 0) {
        fileFilter <- paste(prefix, fileFilter, sep="") 
        if(debug){
          print("prefix added:")
          print(fileFilter)
        }  
      }
      if (!is.na(suffix) && nchar(suffix) > 0) {
        fileFilter <- paste(fileFilter, suffix, sep="") 
        if(debug){
          print("suffix added:")
          print(fileFilter)
        }  
      }
      fileFilter <- paste(fileFilter,"\\.", extension, sep="") 
      if(debug){
        print("fileFilter")
        print(fileFilter)
        print("folder")
        print(folder)
      }  
      
      # Get list of files.
      import.file <- list.files(path = folder, pattern = fileFilter,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = ignore.case, include.dirs = FALSE)
    }
  
  }
    
  if(debug){
    print("import.file")
    print(import.file)
  }  
  
  # Check if files are specified.
  if (any(length(import.file) > 0, !is.na(import.file))) {

    # Autotrim message (function inside loop).
    if(auto.trim){
      message(paste("Auto trim samples:", trim.samples,
                    " invert =", trim.invert))
    }
    
    # Read files.
    for (f in seq(along=import.file)) {

      # Should change to more efficient and simpler 'fread' but
      # problem is that autodetection of colClasses does not always work
      # and it is not possible(?) to set all to character.
      # Use read.table for the time being.
      # Read a file.  
      # tmpdf <- data.table::fread(import.file[f], data.table=FALSE)
      
      # Ensures column names are identical as when read.table was used.
      # Needed since many functions specify columns by name.
      # names(tmpdf) <- make.names(colnames(tmpdf))

      # Read a file.  
      tmpdf <- read.table(import.file[f], header = TRUE,
                          sep = separator, fill = TRUE,
                          na.strings = na.strings,
                          colClasses = "character",
                          stringsAsFactors=FALSE)
      
      # Autotrim datset (message before loop).
      if(auto.trim){
        tmpdf <- trim(data=tmpdf, samples=trim.samples,
                    invert.s=trim.invert, debug=debug)
      }
      
      # Show progress.
      message(paste("Importing (", f, " of ", length(import.file),"): ",
                    import.file[f], sep=""))

      # Check if file path should be saved.
      if(file.name){

        # Add column and save file name.
        tmpdf$File.Name <- basename(import.file[f])
        
      }
      
      # Check if time stamp should be saved.
      if(time.stamp){
        
        # Add column and save file name.
        tmptime <- file.info(import.file[f])
        tmpdf$File.Time <- as.character(tmptime$mtime)
        
      }
      
      # Check if multiple files.
      if(f > 1){
        # Add to result data frame.
        res <- plyr::rbind.fill(res, tmpdf)
      } else {
        # Create result data frame.
        res <- tmpdf
      }
      
    }

    # Autoslim dataset.
    if(auto.slim){
      
      # Autodetect column names to keep fixed.
      fixCol <- colNames(data=res, slim=TRUE, numbered=TRUE, concatenate="|")
      
      # Autodetect column names to stack.
      stackCol <- colNames(data=res, slim=FALSE, numbered=TRUE, concatenate="|")
      
      # Progress.
      message("Auto slim dataset...")
      message(paste("  Stack columns:", stackCol))
      message(paste("  Fix columns:", fixCol))
      
      # Slim require a vector of strings.
      fixCol <- unlist(strsplit(fixCol, "|", fixed = TRUE))
      stackCol <- unlist(strsplit(stackCol, "|", fixed = TRUE))
      
      # Slim data.      
      res <- slim(data=res, fix=fixCol, stack=stackCol,
                  keep.na=slim.na, debug=debug)
      
    }
    
  }

  # Convert common known numeric columns.
  res <- colConvert(data=res)
    
  # Add attributes.
  attr(res, which="import, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
  attr(res, which="import, call") <- match.call()
  attr(res, which="import, date") <- date()
  attr(res, which="import, folder") <- folder
  attr(res, which="import, extension") <- extension
  attr(res, which="import, suffix") <- suffix
  attr(res, which="import, prefix") <- prefix
  attr(res, which="import, import.file") <- import.file
  attr(res, which="import, folder.name") <- folder.name
  attr(res, which="import, file.name") <- time.stamp
  attr(res, which="import, separator") <- separator
  attr(res, which="import, ignore.case") <- ignore.case
  attr(res, which="import, auto.trim") <- auto.trim
  attr(res, which="import, trim.samples") <- trim.samples
  attr(res, which="import, trim.invert") <- trim.invert
  attr(res, which="import, auto.slim") <- auto.slim
  attr(res, which="import, slim.na") <- slim.na
  attr(res, which="import, na.strings") <- na.strings

  return(res)
  
}
