write.dbf <- function(dframe, filename) {

################################################################################
# Function: write.dbf
# Purpose: Write a data frame to the dbf file of an ESRI shapefile
# Programmer: Tom Kincaid
# Date: September 30, 2009
# Revised: July 14, 2014
# Description:
#   This function writes a data frame to a dbf file.
# Arguments:
#   dframe = a data frame to be written to the dbf file
#   filename = name of the dbf file without any extension.
# Results:
#   A data frame composed of either the contents of the single dbf file, when 
#     filename is provided, or the contents of the dbf file(s) in the working
#      directory, when filename is NULL.
# Other Functions Required:
#   writeDbfFile - C function to write a single dbf file or multiple dbf files
################################################################################

# If necessary, strip the file extension from the file name
   if(!is.null(filename)) {
      nc <- nchar(filename)
      if(substr(filename, nc-3, nc) == ".dbf") {
         filename <- substr(filename, 1, nc-4)
      }
   }

# Convert factors to character vectors
   temp <- sapply(dframe, is.factor)
   if(any(temp)) {
      for(i in seq(ncol(dframe))[temp]) {
         dframe[,i] <- as.character(dframe[,i])
         temp <- dframe[,i] == "" | is.na(dframe[,i])
         if(any(temp)) {
            dframe[temp,i] <- " "
         }
      }
      temp <- .Call("writeDbfFile", names(dframe), dframe, filename)
   } else {
      temp <- .Call("writeDbfFile", names(dframe), dframe, filename)
   }

# Return NULL
   invisible(NULL)
}
   