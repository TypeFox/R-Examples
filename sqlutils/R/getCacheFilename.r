#' Returns the complete filepath to the cache file.
#' 
#' @param query the query name.
#' @param dir the directory to save the cache file to.
#' @param ext file extension.
#' @param ... query parameters.
#' @return full filepath to the cached file.
#' @export
getCacheFilename <- function(query, dir=getwd(), ext='csv', ...) {
	parms = getParameters(query)
	parmvals = unlist(list(...))
	filename = paste(dir, '/', query, sep='')
	if(length(parms) > 0) {
		for(i in 1:length(parms)) {
			filename = paste(filename, parms[i], parmvals[parms[i]], sep='.')
		}
	}
	if(nchar(filename) >= 251) {
		warning(paste0('The cached filename is longer than 255 characters. ',
			'This will cause an error on some operating systems. Consider ',
			'specifying your own filename parameter. The filename will be ',
			'truncated to 255 characters.'))
		filename <- substr(filename, 1, 251)
	}
	filename = paste(filename, ext, sep='.')
	return(filename)
}
