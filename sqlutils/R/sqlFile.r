#' Returns the full path to the query or NULL if not found.
#' 
#' @param query the query to find.
#' @return path to the query file.
sqlFile <- function(query) {
	paths <- sqlPaths()
	for(p in paths) {
		f <- paste(p, '/', query, '.sql', sep='')
		if(file.exists(f)) {
			return(f)
		}
	}
	warning(paste(query, ' not found.', sep=''))
	invisible(NULL)
}
