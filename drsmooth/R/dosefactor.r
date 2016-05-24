#' dosefactor
#' @description
#' Returns factorized dose column dose_fac in second column
#' @param dosecolumn   Name of dose column in dataframe.
#' @param data   Input dataframe.
#' @keywords internal

dosefactor <- function (dosecolumn = "", data = NA) {

    data$dose_fac <- as.factor(data[,dosecolumn])
	endcolumn <- ncol(data)
	nexttolastcolumn <- ncol(data)-1
	data <- data[,c(1,endcolumn,2:nexttolastcolumn)]
	return(data)
}

