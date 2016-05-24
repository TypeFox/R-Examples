#' Search paths for SQL repositories.
#' 
#' @param path new path to add. This can a character vector of length greater than 1.
#' @param replace if FALSE, the path(s) will be added to already existing list.
#' @export
sqlPaths <- function(path, replace=TRUE) {
	paths <- unlist(mget("sqlrepos", envir=sqlutils.envir, 
			ifnotfound=list(paste(system.file(package='sqlutils'), '/data', sep=''))))
	if(!missing(path)) {
		path <- normalizePath(path.expand(path), mustWork=FALSE)
		if(replace) {
			paths <- unique(c(path))
		} else {
			paths <- unique(c(path, paths))
		}
		assign("sqlrepos", value=paths, envir=sqlutils.envir)
	}
	return(unname(paths))
}
