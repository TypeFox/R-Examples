read.dbf <- function(filename=NULL) {

################################################################################
# Function: read.dbf
# Purpose: Read the dbf file of an ESRI shapefile
# Programmer: Tom Kincaid
# Date: March 1, 2005
# Last Revised: October 8, 2008
# Description:
#   This function reads either a single dbf file or multiple dbf files.  For 
#   multiple dbf files, all of the dbf files must have the same variable names.
# Arguments:
#   filename = name of the dbf file without any extension.  If filename equals a
#     dbf file name, then that dbf file is read.  If filename equals NULL, then
#     all of the dbf files in the working directory are read.  The default is
#     NULL.
# Results:
#   A data frame composed of either the contents of the single dbf file, when 
#     filename is provided, or the contents of the dbf file(s) in the working
#      directory, when filename is NULL.
# Other Functions Required:
#   readDbfFile - C function to read a single dbf file or multiple dbf files
################################################################################

# If necessary, strip the file extension from the file name

   if(!is.null(filename)) {
      nc <- nchar(filename)
      if(substr(filename, nc-3, nc) == ".dbf") {
         filename <- substr(filename, 1, nc-4)
      }
   }

# Read the dbf file

   dbffile <- .Call("readDbfFile", filename)
   if(is.null(dbffile[[1]]))
      stop("\nAn error occurred while reading the dbf file(s) in the working directory.")

# Convert character vectors to factors

   ind <- sapply(dbffile, is.character)
   if(any(ind)) {
      for(i in (1:length(dbffile))[ind])
         dbffile[,i] <- as.factor(dbffile[,i])
   }

# Assign class to the output data frame

   attr(dbffile, "class") <- c("SurveyFrame", "data.frame")

# Return the data frame

   dbffile
}
   