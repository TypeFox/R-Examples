"getsonde" <- function(filename, datakey="------", varkey=" Time", unitkey="  sec")
{
#
# The file structure is that there are a bunch of header lines
# followed by some delimiting string that indicates the NEXT line
# is the start of the data ...
#

#
# Copyright 2001,2002 Tim Hoar, Eric Gilleland, and Doug Nychka
#
# This file is part of the RadioSonde library for R and related languages.
#
# RadioSonde is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# RadioSonde is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RadioSonde; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

# must parse the file initially and find:
#   a) the line with the column headers 
#   b) the line with the units for each column
#   c) the line that precedes the data 

	foo    <- scan(filename,what="character",sep="\r",quiet=TRUE)
        nlines <- length(foo)
        inds   <- seq(1:nlines)

#
# Decode the column names / variable fields, whatever
#

	if ( ! is.null(varkey) ) {
	   nameindex <- grep(paste("^",varkey,sep=""),foo)
           if ( length(nameindex) == 0 ){
           stop(paste("(getsonde): Variable-name String <",varkey,"> not found.",sep=""))
           }   
	} else {
           stop("(getsonde): Must be a non-NULL 'varkey'")
	}

        if ( length(nameindex) != 1 ) {
           stop("(getsonde): could not find a unique match for the variable-name string")
        }

	col.names <- tolower(scan(filename, what="", skip=nameindex-1,
                             nlines=1, strip.white=TRUE, quiet=TRUE))

#
# Figure out how many lines to skip and read into a dataframe
#

	if ( ! is.null(datakey) ) {
	   dataindex <- grep(paste("^",datakey,sep=""),foo)
           if ( length(dataindex) == 0 ){
           stop(paste("(getsonde): Data String <",datakey,"> not found.",sep=""))
           }   
	} else {
           stop("(getsonde): Must be a non-NULL 'datakey'")
	}

        if ( length(dataindex) != 1 ) {
           stop("(getsonde): could not find a unique match for the data string")
        }

        dataframe <- read.table(filename,skip=dataindex,col.names=col.names)

	missingones <- (dataframe$press == 9999. | dataframe$rh == 999. | dataframe$dewpt == 999.)

        dataframe[missingones,] <- NA

#
# If there are units, try to decode them
#

	if ( ! is.null(unitkey) ) {
	   unitindex <- grep(paste("^",unitkey,sep=""),foo)
           if ( length(unitindex) == 0 ){
           stop(paste("(getsonde): Unit String <",unitkey,"> not found.",sep=""))
           }   

	   units <- scan(filename, what="", skip=unitindex-1, nlines=1, strip.white=TRUE, quiet=TRUE)
           attr(dataframe,"units") <- units

	}

#
# Store the header information as ancillary data
#

        attr(dataframe,"metadata") <- foo[1:dataindex]

	return(dataframe)
}
