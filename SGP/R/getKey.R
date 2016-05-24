`getKey` <- 
function(sgp_object) {
	if (is.SGP(sgp_object)) {
		if ("YEAR_WITHIN" %in% names(sgp_object@Data)) return(c("VALID_CASE", "CONTENT_AREA", "YEAR", "ID", "YEAR_WITHIN")) else return(c("VALID_CASE", "CONTENT_AREA", "YEAR", "ID"))
	}
	if (is.data.table(sgp_object)) {
		if ("YEAR_WITHIN" %in% names(sgp_object)) return(c("VALID_CASE", "CONTENT_AREA", "YEAR", "ID", "YEAR_WITHIN")) else return(c("VALID_CASE", "CONTENT_AREA", "YEAR", "ID"))
	}
} ### END getKey
