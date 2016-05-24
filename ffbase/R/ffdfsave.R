#' Save a ffdf data.frame in directory
#'
#' \code{ffdfsave} saves a ffdf data.frame in the given filename (.rdata) and stores all
#' \code{ff} columns in a subdirectory with the name "<filename>_ff". Each column
#' will be named "<columnname>.ff".
#' A saved ffdf data.frame is a .rdata file and can be loaded with the \code{load} function
#' Deprecated, the preferred method is \code{\link{save.ffdf}}
#' @rdname pkg-deprecated
#' @export
#' @param dat \code{ffdf} data.frame, to be saved
#' @param filename path where .rdata file will be save and <filename>_ff directory will be created
ffdfsave <- function(dat, filename){
   .Deprecated("Use save.ffdf")
   
   datname <- deparse(substitute(dat))
   
   # create a sub directory with "<filename>_ff"
   dirnm <- sub("(\\..+)?$","_ff",filename)
   dir.create(dirnm, showWarnings=FALSE)
   
   for (colname in names(dat)){
      ffcolname <- paste(file.path(dirnm, colname), "ff", sep=".")
      ffcol <- dat[[colname]]
      # set filename of ff file to '<columnname>.ff' in sub directory
      filename(ffcol) <- ffcolname      
   }
   
   # close all ff files...
   close(dat)
   assign(datname, dat)
   # save ffdf with original name to <filename>
   save(list=datname, file=filename)
}