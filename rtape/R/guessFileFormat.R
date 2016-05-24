#' Automatically pick proper tape file format.
#'
#' This function guesses the tape file format from the file header or assumes default (gzip) if given non-existing file name. Main package functions use this routine to automatically setup file format; if you really need to control it, see \code{\link{makeFileFormat}}.
#'
#' @param fName Name of the the file to guess format of; if the file is not-existing, the function returns default file format.
#' @return The function to be passed to the \code{fileFormat*} arguments of other \code{rtape}'s functions.
#' @author Miron B. Kursa \email{M.Kursa@@icm.edu.pl}

guessFileFormat<-function(fName){
 stopifnot(length(fName)==1)
 if(!file.exists(fName)) return(makeFileFormat())
 readBin(fName,what="int",n=6,size=1,signed=FALSE)->header
 if(identical(header[1:2],c(0x1fL,0x8bL))) return(makeFileFormat("gz"))
 if(identical(header[1:6],c(0xfdL,0x37L,0x7aL,0x58L,0x5aL,0x00L))) return(makeFileFormat("xz"))
 if(identical(header[1:3],c(0x42L,0x5aL,0x68L))) return(makeFileFormat("bz"))
 stop("Tape is corrupted!")
}
