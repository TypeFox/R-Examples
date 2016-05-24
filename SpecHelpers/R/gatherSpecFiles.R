#' Process LoggerPro Spectral Files into a Data Frame
#' 
#' This function will go through all the files of a specified format in a
#' directory and convert them into a data frame with one column containing the
#' wavelength information and the other columns the absorbances of each sample
#' (file).  The file names are used to create the column names in the data
#' frame.  Optionally, non-integer wavelengths in the file can be combined to
#' give integer wavelengths. Keep in mind that this function specifically
#' modifies formats written by LoggerPro.  Each format, as it comes from
#' LoggerPro, has various amounts of crap in it which has to be removed or
#' modified.
#' 
#' All files of a given extension in the directory will be processed, so make
#' certain there are no extra files in the directory.  The files will be
#' modified and written back out as .csv files so look for the number of files
#' in the directory to double.  In the case of csv files, the original csv
#' files will be overwritten.  These files have no header row.
#' 
#' @param type A character string giving the type of files to be processed.
#' Currently, either "txt", "csv" or "cmbl" extensions can be processed.
#'
#' @param intLambda Logical.  If TRUE, non-integer wavelengths that round to
#' the same value will be combined and averaged and reported as integer values.
#'
#' @param \dots Other parameters to be passed downstream.  Currently none
#' possible.
#' 
#' @return A data frame containing the wavelengths in the first column and the
#' absorbances in the other columns, one column per file, with column names
#' generated from the file names.
#'
#' @author Bryan A. Hanson, DePauw University
#' @keywords utilities
#' @export
#' @importFrom utils read.csv
#'
gatherSpecFiles <-
function(type = "txt", intLambda = FALSE, ...) {

	message("Are you working on a copy?")
	answer <- substr(readline("Continue (y/n)?  "), 1L, 1L)
		if (answer == "n" | answer == "N") {
			cat("cancelled by user\n")
			return(invisible())
			}
	if (answer == "y" | answer == "Y") cat("Onward then!\n")

	if (type == "csv") message("Files are being overwritten...")
	
	if (type == "cmbl") {
		files <- list.files(pattern = "\\.(cmbl|CMBL)")
		files.noext <- substr(basename(files), 1, nchar(basename(files)) - 4)
		out.files <- paste(files.noext, "csv", sep = "")
	
		for (i in 1:length(files)) {
			cmbl2csv(in.file = files[i], out.file = out.files[i])
			}
		}
		
	if (type == "txt") {
		files <- list.files(pattern = "\\.(txt|TXT)")
		files.noext <- substr(basename(files), 1, nchar(basename(files)) - 3)
		out.files <- paste(files.noext, "csv", sep = "")
	
		for (i in 1:length(files)) {
			txt2csv(in.file = files[i], out.file = out.files[i])
			}
		}

	if (type == "csv") { # read 'em and remove header
		files <- list.files(pattern = "\\.(csv|CSV)")
	
		for (i in 1:length(files)) {
			df <- read.csv(files[i])
			write.table(df, file = files[i], row.names = FALSE,
				col.names = FALSE, quote = FALSE, sep = ",")
			}
		}

	if (type == "sstab") {
		files <- list.files(pattern = "\\.(txt|TXT)")
		files.noext <- substr(basename(files), 1, nchar(basename(files)) - 3)
		out.files <- paste(files.noext, "csv", sep = "")
	
		for (i in 1:length(files)) {
			sstab2csv(in.file = files[i], out.file = out.files[i])
			}
		}
	
	if (intLambda) avgLambda() # aggregate wavelengths into whole numbers and mean abs
	
	df <- gatherCsv(...) # combine the many into one
	
	}

