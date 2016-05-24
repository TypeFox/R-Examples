#' R functions, a class, and methods for working with SaTScan stand-alone software.
#' 
#' rsatscan provides a suite of functions that allows you to easily write SaTScan parameter
#' files in the OS, run SaTScan in the OS, and read the output files that SaTScan generates.
#' 
#' The parameter files are constructed in R using the \code{ss.options} function and written 
#' to the OS using the \code{write.ss.prm} function.  SaTScan is run using the \code{satscan}
#' function.  The \code{satscan} function returns a  \code{satscan-class} object that has a
#' slot for every possible file that SaTScan makes, plus one for the parameter file you used
#' to generate the output.
#' 
#' The package also includes  \code{write.???} functions which will write case, control, 
#' geography, population, etc., files in the format expected by SaTScan, if you happen to have
#' them (or make them) in R and want to write them into the OS for SaTScan to use.
#' 
#' There are \code{summary} and \code{print} methods for \code{satscan-class} 
#' objects.  There are also \code{plot} methods in the \code{sp} package, which
#' can be used if the \code{rgdal} package and \code{sp} packages are
#' installed and SaTScan generated a shapefile.
#' 
#' Currently the package works with SaTScan >= 9.2 and has been tested on Windows 7 and 
#' Ubuntu 14.04.1.  Please contact the author if you find success or trouble on other OSes. 
#' 
#' @docType package
#' @name rsatscan
NULL