################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 29.08.2015: Added importFrom.
# 15.12.2014: Changed parameter names to format: lower.case
# 30.09.2014: Check if package is loaded to avoid error in path.package.
# 01.10.2013: First version.

#' @title Get Allele Frequency Database
#'
#' @description
#' Gives access to allele frequency databases.
#'
#' @details
#' The function provides access to allele frequency databases stored in
#' the file database.txt in the package directory.
#' It returns the specified allele frequency database.
#' 
#' @param db.name.or.index string or integer specifying the database.
#' If NULL a vector of available databases is returned.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame with allele frequency database information.
#' If no matching database or database index is found NA is returned.
#' If the database was not found NULL is returned.
#' 
#' @keywords internal
#' 
#' @export
#'  
#' @importFrom utils read.delim
#' 
#' @examples
#' # Show available allele frequency databases.
#' getDb()

getDb <- function(db.name.or.index=NULL, debug=FALSE) {

  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }

  .separator <- .Platform$file.sep # Platform dependent path separator.
  
  # LOAD DATABASE INFO  #######################################################

  # Get package path.
  if(require(strvalidator)){
    packagePath <- path.package("strvalidator", quiet = FALSE)
  } else {
    warning("Package path for strvalidator not found!")
    return(NULL)
  }
  subFolder <- "extdata"
  fileName <- "database.txt"

  # Create complete file path.
  filePath <- paste(packagePath, subFolder, fileName, sep=.separator)
  
  .db <- read.delim(file=filePath, header = TRUE, sep = "\t", quote = "\"",
                         dec = ".", fill = TRUE, stringsAsFactors=FALSE)
  
  # Available databases. Must match else if construct.
  databases<-unique(.db$Database)
  
	# Check if NULL
	if (is.null(db.name.or.index)) {

		db <- databases

	# String provided.
	} else {

		# Check if number or string.
		if (is.numeric(db.name.or.index)) {

			# Set index to number.
			index <- db.name.or.index

		} else {

			# Find matching database index (case insensitive)
			index <- match(toupper(db.name.or.index),toupper(databases))

		}

		# No matching database.
		if (is.na(index)) {
			
			db <- NA

		# Assign matching database information.
		} else {
		  
		  db <- .db[.db$Database==databases[index], ]
      
		} 

	}
  
  return(db)

}
