objSearch <- function(sep = "\t", path = NULL, compare = TRUE)
{
    Search <- search()
	if (isTRUE(compare)) {
		oldSearch <- getTemp(".guiObjSearchCache", default = "")
		## Compare both versions
		if (length(Search) != length(oldSearch) || !all(Search == oldSearch)) {
			## Keep a copy of the last version in SciViews:TempEnv
			assignTemp(".guiObjSearchCache", Search)
			Changed <- TRUE
		} else Changed <- FALSE
	} else Changed <- TRUE
    if (is.null(path)) {  # Return result, as a single character string with sep
		if (Changed) {
			if (!is.null(sep)) Search <- paste(Search, collapse = sep)
			return(Search)
		} else return("")
	} else {  # Write to a file called 'Search.txt' in this path
		file <- file.path(path, "Search.txt")
		if (Changed) {
			if (is.null(sep)) sep <- "\n"
			cat(Search, sep = sep, file = file)
		}
		return(invisible(Changed))
	}
}
