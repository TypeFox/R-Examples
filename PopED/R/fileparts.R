#' MATLAB fileparts function
#' 
#' Get the various parts of a file with path string.
#' 
#' @param filename.with.path A string of a filename with a path
#' 
#' @return A list with the following components:
#' \item{pathname}{The path name}
#' \item{filename}{The file name}
#' \item{fileext}{The file extension}
#' 
#' @note This is a modified version of the same function in \code{\link[matlab]{fileparts}}
#' 
#' @family MATLAB
#' @example tests/testthat/examples_fcn_doc/examples_fileparts.R
#' @export
#' @keywords internal
## Function written to match MATLAB function
## Author: Andrew Hooker

fileparts <- function(filename.with.path){
    pathname <- dirname(filename.with.path)
    filename <- basename(filename.with.path)
    fileext <- gsub(".*(\\.[^\\.]*)$","\\1",filename)
    filename <- gsub("(.*)(\\.[^\\.]*)$","\\1",filename)
    if(fileext==filename) fileext <- ""
    return(list(pathname=pathname,filename=filename,fileext=fileext))
}
