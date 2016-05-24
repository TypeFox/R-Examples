#' @title Case Attributes from Xes-file
#' @description Extracts case attributes from a xes-file.
#' @param xesfile Reference to a .xes file, conforming to the xes-standard.
#' @seealso \url{http://www.xes-standard.org/}
#' @export case_attributes_from_xes

case_attributes_from_xes <- function(xesfile = file.choose()) {
	parsed_xes <- parseXES(xesfile)
	n_case_att <- length(parsed_xes$trace.att)
	n_event_att <- length(parsed_xes$event.att)
	n_cases <- length(parsed_xes$traces$'concept:name')
	caseids <- parsed_xes$traces$'concept:name'
	result <- data.frame(stringsAsFactors = FALSE)
	data <- parsed_xes$traces
	for(i in 1:length(parsed_xes$traces[[1]])){
		for(j in 1:n_case_att){
			if(i > length(data[[j]]))
				result[i,j] <- NA
			else
				result[i,j] <- data[[j]][[i]]
		}
	}
	for(i in 1:n_case_att)
		colnames(result)[i] <- paste("case",gsub(":",".", names(parsed_xes$traces)[i]), sep = "_")
	result <- tbl_df(result)
	return(result)
}
