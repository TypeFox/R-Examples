#' Non-existent files from a vector of filenames
#' 
#' This function receives a vector of filenames as parameter, 
#' and returns a vector of non-existent files among them.
#' 
#' @param filenames character A vector of filenames
#' @return A character vector of file paths that do not exist.
#' @examples 
#' \dontrun{
#' nonExistentFiles(R.home())
#' nonExistentFiles(sapply(1:5, function(i) tempfile()))
#' nonExistentFiles(sapply(1:5, function(i) tempdir()))
#' nonExistentFiles(c("/tmp/f3412lds43289ajkfdlsa", R.home())) == "/tmp/f3412lds43289ajkfdlsa"
#' }
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
nonExistentFiles = function(filenames) {
	filenames[!file.exists(filenames)]
}


#' Stop when any file does not exist 
#' 
#' @param files character vector. File paths you want to check.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
checkFileExist = function(files) {
	miss_files = nonExistentFiles(files)
	if(length(miss_files) > 0) {
		stop(paste("Non-existent files:", paste(miss_files, collapse=", ")))
	}
}
