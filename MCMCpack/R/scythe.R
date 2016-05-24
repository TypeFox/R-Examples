##########################################################################
## Scythe Inter-Operation Functions
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


# writes a matrix out to an ASCII file that can be read by Scythe.
# it puts the number of rows and columns in the first row
# followed by the data.
#
# ADM 1/29/2003

"write.Scythe" <-
  function(outmatrix, outfile = NA, overwrite=FALSE) {
    outmatrix <- as.matrix(outmatrix)
   
    if(is.na(outfile)) {
      stop("Please specify a file name in the write.Scythe() call.\n")
    }
    if(overwrite==FALSE & file.exists(outfile)) {
      cat("File already exists in the write.Scythe() call.\n")
      stop("Either delete the file, or flip the overwrite switch.\n")
    }
   
    outfile <- file(outfile, "w")
    cat(dim(outmatrix), "\n", file=outfile)
    write.table(outmatrix, file=outfile,
                row.names=FALSE, col.names=FALSE, quote=FALSE)
    close(outfile)
    return(0)
  } 


# reads in a matrix from an ASCII file written by Scythe.
# the number of rows and columns should be in the first row followed
# by the data.
#
# Kevin Rompala 5/1/2003
# fixed by ADM 7/25/2004

"read.Scythe" <-
  function(infile = NA) {
    
    if(is.na(infile)) {
      stop("Please specify a file name in the read.Scythe() call.\n")
    }
    if(!file.exists(infile)) {
      stop("Specified source file does not exist in read.Scythe() call.\n")
    }

    infile <- file(infile, "r")
    dimensions <- scan(file=infile,n=2)
    inputdata <- scan(file=infile)
    close(infile)
    hold <- matrix(data=inputdata,
           nrow=dimensions[1], ncol=dimensions[2], byrow=TRUE)
    return(hold) 
  }
