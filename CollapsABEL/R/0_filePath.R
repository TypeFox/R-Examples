#' An S4 class to represent a file path
#' 
#' This class comes with a validation function, 
#' making sure that the file exists.
#' 
#' @slot path character, file or dir path
#' 
#' @author Kaiyin Zhong, Fan Liu
.FilePath = setClass("FilePath", representation(path = "character"), 
		validity = function(object) {
			missing_files = nonExistentFiles(object@path)
			if(length(missing_files) > 0) {
				paste("Missing file", missing_files)
			} else {
				TRUE
			}
		})

#' Constructor for FilePath class 
#' 
#' @param s character, path to file or dir
#' @return FilePath object
#' @examples 
#' \dontrun{
#' fp = filePath(R.home())
#' dirName(fp) == dirname(fp@path)
#' baseName(fp) == basename(fp@path)
#' }
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
filePath = function(s) {
	new("FilePath", path = s)
}

#' Directory name of a file path
#'  
#' @name dirName
#' 
#' @param fp FilePath object
#' @return character vector of directories 
#' @examples 
#' \dontrun{
#' fp = filePath(R.home())
#' dirName(fp)
#' }
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
setGeneric("dirName",
		function(fp) {
			standardGeneric("dirName")
		})

#' @rdname dirName
#' @docType methods
#' @export
setMethod("dirName",
		signature(fp = "FilePath"),
		function(fp) {
			dirname(fp@path)
		})

#' Basename of a FilePath object
#' 
#' @name baseName
#' 
#' @param fp character, file paths.
#' @return character vector of basenames
#' @examples
#' \dontrun{
#' fp = filePath(R.home())
#' baseName(fp)
#' }
#' 
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
setGeneric("baseName",
		function(fp) {
			standardGeneric("baseName")
		})


#' @rdname baseName 
#' @export 
setMethod("baseName",
		signature(fp = "FilePath"),
		function(fp) {
			basename(fp@path)
		})



