#' Combine csv Files Containing Spectral Data into a Data Frame
#' 
#' This function processes csv files containing two columns, wavelength and
#' absorbance (or intensity etc), into a data frame, which is then written out as a csv file.  The
#' files should have no header row.
#' 
#' It is assumed that the csv files have already been cleaned up so that they
#' contain only wavelength and absorbance data.  The wavelength data column
#' must be the same in all the files (as they would be if they came from the
#' same instrument with the same settings).
#' 
#' @return A data frame containing the wavelengths in the first column and the
#' absorbances in the other columns, one per file, with the file name
#' generating the column name.  The data frame is written out in a file called
#' "All Spec Files.csv".
#'
#' @author Bryan A. Hanson, DePauw University
#' @seealso \code{\link{gatherSpecFiles}} which is the function the user should
#' call.
#' @keywords utilities
#' @export
#' @importFrom utils read.csv write.csv
gatherCsv <-
function() {

	files <- list.files(pattern = ".csv|.CSV")
	names <- substr(basename(files), 1, nchar(basename(files)) - 4)
	wave <- read.csv(files[1], header = FALSE)
	df <- wave

	for (i in 2:length(files)) {
		temp <- read.csv(files[i], header = FALSE)
		temp <- temp[,2]
		df <- data.frame(df, temp)
		}

	colnames(df) <- c("Wavelength", names)
	write.csv(df, "All Specs.csv", row.names = FALSE)
	df
	}

