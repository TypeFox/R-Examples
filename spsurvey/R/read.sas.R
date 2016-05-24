read.sas <- function(filename, libname=NULL, xport=FALSE,
   sascmd="/Program Files/SAS/SAS 9.1/sas.exe") {

################################################################################
# Function: read.sas
# Purpose: Read SAS datasets or a SAS XPORT (transport) file
# Programmer: Tom Kincaid
# Date: May 9, 2007
# Revised: March 9, 2010
# Description:
#   This function reads either a SAS dataset or a SAS XPORT (transport) file and
#   creates a data frame.
# Arguments:
#   filename = if xport equals TRUE, a character string giving the full path to
#     the SAS XPORT file, which must include the file extension.  If xport
#     equals FALSE, either a character string giving the the name of a dataset
#     in the SAS library or a vector of character strings giving the names of
#     datasets in the SAS library, where the dataset names cannot exceed eight
#     characters in length and do not include the file extension.
#   libname = a character string defining the SAS library, which is usually a
#     directory reference.   If xport equals FALSE and the dataset(s) named in
#     argument filename do not reside in the working directory, then this
#     argument is required.  The default value is NULL.
#   xport = a logical value indicating whether the input file is a SAS XPORT
#     file.  The default value is FALSE.
#   sascmd = a character string giving the full path to SAS executable.  This
#     argument is required only when xport equals FALSE.  The default value is
#     "C:/Program Files/SAS/SAS 9.1/sas.exe".
# Results:
#   Either a single data frame or a list of data frames.
# Other Functions Required:
#   read.ssd - function in the foreign library that reads a SAS dataset and
#     creates a data frame
#   read.xport - function in the foreign library that reads a SAS XPORT file and
#     creates a data frame
# Example:
#   MySasFile <- read.sas("mysasfil", 
#                         "C:/Documents and Settings/auser/My Project")
################################################################################

# Read the SAS XPORT file

   if(xport) {
      df <- read.xport(filename)

# Read the SAS dataset(s)

   } else if(is.null(libname)) {
      df <- read.ssd(getwd(), filename, sascmd=sascmd)

   } else {
      df <- read.ssd(libname, filename, sascmd=sascmd)
   }

# Return the data frame(s)

   return(df)
}
