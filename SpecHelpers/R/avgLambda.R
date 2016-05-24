#' Convert Wavelengths to Integer Values and Average Corresponding Absorbances
#' 
#' This function reads a csv file containing columns of wavelength and
#' absorbances.  It rounds the wavelengths to integers and replaces the
#' absorbances corresponding to the rounded values with their averages.  NOTE
#' THAT ALL THE csv FILES IN THE CURRENT DIRECTORY ARE PROCESSED AND REPLACED
#' WITH THE MODIFIED FILEs.  You should use this function on a copy of the
#' directory.
#' 
#' 
#' @return The original files are overwritten with the modified files.
#'
#' @author Bryan A. Hanson, DePauw University
#' @seealso \code{\link{gatherSpecFiles}} which is the function the user should
#' call.
#' @keywords utilities
#' @export
#' @importFrom utils read.csv write.table
#' @importFrom stats aggregate na.pass
avgLambda <-
function() {

	files <- list.files(pattern = ".csv|.CSV")
	for (i in 1:length(files)) {
		temp <- read.csv(files[i], header = FALSE)
		#cat("Ready to round wavelengths ", files[i], "\n")
		temp[,1] <- round(temp[,1])
		#cat("Ready to average ", files[i], "\n")
		temp <- aggregate(data = temp, V2~V1, mean, na.action = na.pass) # what a time saver!
		write.table(temp, file = files[i], row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
		}
	}

