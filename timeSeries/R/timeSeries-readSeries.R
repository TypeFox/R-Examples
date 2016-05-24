#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  ../../COPYING

################################################################################
# FUNCTION:                 DESCRIPTION:
#  readSeries                Reads a CSV file and creates a 'timeSeries'
################################################################################


# DW:
# I think we should add a similar function for writeSeries() using 
# write.table(). Proceed in the same way as in the case of the read
# function.


# ------------------------------------------------------------------------------


readSeries <-
    function(file, header = TRUE, sep = ";", zone = "", FinCenter = "",
    format, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Reads from a spreadsheet and creates a 'timeSeries' object

    # Arguments:
    #   file - the name of the file which the data are to be read 
    #       from. Each row of the table appears as one line of the 
    #       file. If it does not contain an absolute path, the file 
    #       name is relative to the current working directory, 
    #       getwd(). Tilde-expansion is performed where supported. 
    #       As from R 2.10.0 this can be a compressed file.
    #   header - a logical value indicating whether the file contains 
    #       the names of the variables as its first line. If missing, 
    #       the value is determined from the file format: header is 
    #       set to TRUE if and only if the first row contains one fewer 
    #       field than the number of columns.
    #   sep - he field separator character. Values on each line of 
    #       the file are separated by this character. If sep = "" (the 
    #       default for read.table) the separator is ?white space?, 
    #       that is one or more spaces, tabs, newlines or carriage 
    #       returns.
    #   zone - the time zone or financial center where the data were
    #       recorded.
    #   FinCenter - a character with the the location of the
    #       financial center named as "continent/city". By default
    #       an empty string which means that internally "GMT" will
    #       be used.
    #   format - the format of the timestamps as recoreded in the
    #       first column of the data in the..
    #   ... - optional arguments passed to the function read.table().
    
    # Value:
    #   Returns a S4 object of class 'timeSeries'.

    # Notes:
    #   Note we expect that the header of the spreadsheet file in
    #   the first cell holds the time/date format specification!

    # FUNCTION:

    # Read Data:
    df <- read.table(file = file, header = header, sep = sep,
        check.names = FALSE, ...)

    # Get 'timeDate' from first column with header specifying the format
    charvec <- as.character(df[[1]])
    if (missing(format)) format <- names(df)[1]
    td <- try(timeDate(charvec = charvec, format = format, zone = zone,
        FinCenter = FinCenter), silent=TRUE)
      
    # DW: 2014-09-16
    # If sep=";" fails try with sep=",":
    if (sep ==";" && class(td) == "try-error") {
        return(readSeries(file, header = header, sep = ",", zone = zone, 
          FinCenter = FinCenter, ...))
    }

    # If format provided in file or with format argument, try to guess it
    if (all(is.na(td)))
        warning("Conversion of timestamps to timeDate objects produced only NAs.
  Are you sure you provided the proper format with argument 'format'
  or in the header of your file ?")

    # Extract data
    data <- as.matrix(df[-1])

    # Create Time Series from Data Frame:
    ans <- timeSeries(data = data, charvec = td)

    # Return Value:
    ans
}


################################################################################

