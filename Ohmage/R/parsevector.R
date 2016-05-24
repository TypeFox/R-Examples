parsevector <- function(x){
	#the main purpose of this function is to map JSON null values to NA.
	ifelse(is.null(x) || x == "NOT_DISPLAYED" || x == "SKIPPED", return(NA), return(x))
}
