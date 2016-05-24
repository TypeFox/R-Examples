#' Utility Functions to Clean and Convert Spectral Files to csv
#' 
#' These functions clean out extraneous information from exported spectral data
#' files and then write them out in csv format.  \code{txt2csv} and \code{cmbl2csv} handle
#' files exported by LoggerPro software.  \code{sstab2csv} handles files exported by
#' Spectra Suite software.  Not directly called by the user.
#' 
#' Extraneous text at the beginning of the file is removed.  In the case of
#' cmbl files, lines containing "Z2" or ">" are removed.  Absorbances marked as
#' "Z1" are replaced with zero.  The data are initially in one long column; the
#' wavelength and absorbances are reunited into two columns.
#' 
#' @aliases txt2csv cmbl2csv sstab2csv
#'
#' @param in.file The name of the input file.
#'
#' @param out.file The name of the output file.
#'
#' @return A modifed file in csv format.
#'
#' @author Bryan A. Hanson, DePauw University.
#'
#' @seealso \code{\link{gatherSpecFiles}} which is the function the user should
#' call.
#' @keywords utilities
#' @importFrom utils write.table read.table
#' @export
#' 
txt2csv <-
function(in.file = "", out.file = "") {
	
	df <- read.table(in.file, skip = 7, sep = "\t")
	write.table(df, file = out.file, row.names = FALSE,
		col.names = FALSE, quote = FALSE, sep = ",")

	}

