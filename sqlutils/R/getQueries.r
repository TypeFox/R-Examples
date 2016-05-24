#' Returns a list of available queries in the current repository.
#' 
#' @export
getQueries <- function() {
	paths <- sqlPaths()
	files <- character()
	for(p in paths) {
		files = c(files, list.files(path=p, pattern="*.sql"))
	}
	return( substr(files, 0, nchar(files)-4) )
}
