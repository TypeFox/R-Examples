#' Function for working with cached queries.
#' 
#' This will first look in the given directory for a CSV or Rda version of the file, if
#' it exists, that will be read and returned. Otherwise it will execute the query
#' and then saves a CSV or Rda file.
#' 
#' @param dir the directory to save and load cached data files. Defaults to the
#'        current working directory (i.e. \code{\link{getwd}}.
#' @param filename the filename of the cached data file.
#' @param query the query to execute.
#' @param maxLevels the maximum number of levels a factor can have before being
#'        converted to a character vector.
#' @param ... other parameters passed to the \code{\link{execQuery}} function including
#'        query parameters.
#' @param format either \code{csv} for comma separated value files or \code{rda} for R data files.
#' @return a data frame.
#' @export
cacheQuery <- function(query=NULL, dir=getwd(), 
					   filename=getCacheFilename(query=query, dir=dir, ext=format, ...), 
					   format='rda', 
					   maxLevels=20, 
					   ...) {
	if(file.exists(filename)) {
		message(paste("Reading from cached query file: ", filename, sep=''))
		if(tolower(format) == 'rda') {
			load(filename)
		} else if(tolower(format) == 'csv') {
			df = read.csv(filename)			
		} else {
			stop('Unsupported format type.')
		}
		df = recodeColumns(df, maxLevels)
	} else {
		message(paste("Executing ", query, " and saving to ", filename, sep=''))
		df = execQuery(query=query, maxLevels=maxLevels, ...)
		if(tolower(format) == 'rda') {
			save(df, file=filename)
		} else if(tolower(format) == 'csv') {
			write.csv(df, filename, row.names=FALSE)
		} else {
			stop('Unsupported format type.')
		}
	}
	return(df)
}
