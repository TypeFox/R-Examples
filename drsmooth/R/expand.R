#' @name expand
#' @aliases expand
#' @title Expand Summarized Dichotomous Data for drsmooth
#' @usage
#' expand(dosecolumn = "", 
#'       targetcolumn = "", 
#'       ncolumn = "",
#'       data = NA,
#'       outputfilename = "")
#' @description
#' This function expands summarized dichotomous data for input into drsmooth.
#' @details
#' The input data frame to this function is expected to contain a numerical dose column,
#' identified as 'dosecolumn,' a numerical outcome count column, identified as 'targetcolumn,'
#' and a numerical N column with the total number of cases tested in the dose group,
#' identified as 'ncolumn.'  See drsmooth included data DIdata for an example.
#' This function takes such a summarized dataframe and returns an expanded case-by-case replicate df,
#' which is suitable as an input to the drsmooth function, with data_type = "dichotomous".
#' If an outputfilename is provided, the resulting df is also written to the working directory.
#' @param dosecolumn  Character string, name of dose column to be tested.
#' @param targetcolumn  Character string, name of response column to be tested.
#' @param ncolumn  Character string, name of N column of total cases per dose.
#' @param data  Input dataframe.
#' @param outputfilename  Character string, name of .csv file to be written out.
#' @return
#' A dataframe is returned containing one row for each case in each dose.
#' If an outputfilename is provided, the df is also written out in a .csv file.
#' An expanded dataframe or file such as this one is appropriate as
#' input to the drsmooth function when data_type = "dichotomous".
#' @examples
#' # Generates an expanded dataframe from the included data DIdata:
#' data(DIdata)
#' DIdata_expanded <- expand(dosecolumn = "Dose", targetcolumn = "Tumor", ncolumn = "n", data = DIdata)
#' @export

# AB 3_10_15 added all necessary documentation for packaging, changed target variable
# and data names, and parameter order to be consistent with the rest of the code in drsmooth.
# Removed data duplication in 'x' object. Made written output file an option, otherwise
# return the expanded dataframe to the current environment for calling drsmooth.
# Removed hard-coding of column headings; retained original column headings with minimum
# impact to existing piecewise dataframe build.

expand <- function (dosecolumn="", targetcolumn="", ncolumn="", data=NA, outputfilename="") {
    
	# x <- data
	dosevariable <- data[,dosecolumn]
	nvariable <- data[,ncolumn]
	responsevariable <- data[,targetcolumn]

	dose <- rep(dosevariable[1],nvariable[1])
	response <- c(rep(0,nvariable[1]-responsevariable[1]),rep(1,responsevariable[1]))

  # Retain original column names
  currentdata <- data.frame(dose, response)
  names(currentdata)[1] <- dosecolumn
  names(currentdata)[2] <- targetcolumn
  
	for (i in 2:nrow(data)) {
		dose <- rep(dosevariable[i],nvariable[i])
		response <- c(rep(0,nvariable[i]-responsevariable[i]),rep(1,responsevariable[i]))
		adddata <- data.frame(dose,response)
		names(adddata)[1] <- dosecolumn
		names(adddata)[2] <- targetcolumn
		currentdata <- rbind(currentdata,adddata)
		}
  
  return(currentdata)
	if (outputfilename != "") {utils::write.csv(currentdata, outputfilename, row.names=FALSE)}

}
