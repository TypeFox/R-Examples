#' Returns the parameters that must be set for the given query.
#' 
#' @param query the query name.
#' @return list of parameter names.
#' @export
getParameters <- function(query) {
	sql = getSQLRaw(query)
	pos = gregexpr(":", sql)
	results = character()
	if(pos[[1]][1] > 0) {
		for(i in seq(1, length(pos[[1]]), by=2)) {
			results = c(results, (substr(sql, pos[[1]][i]+1, pos[[1]][i+1]-1)) )
		}
	}
	return(unique(results))
}
