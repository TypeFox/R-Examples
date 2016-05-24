#' Extract the records from an EpiJson Object
#' 
#' Get the records as a list from an EpiJSON Object
#' @param x An EpiJSON object (class ejObject)
#' @note This is functionally equivelent to calling x$records (which may be 
#'  quicker as therer is no object type checking) but open to change in later 
#'  versions.
#' @return A list of records
#' @export 
ejRecords <- function(x){
	if(class(x) != "ejObject")
		stop("The object supplied to records must be an ejObject")
	x$records
}

#' Extract the attributes from an EpiJson Object
#' 
#' Get the attributes as a list from an EpiJSON Object
#' @param x An EpiJSON object (class ejObject, ejRecord or ejEvent)
#' @return A list of attributes
#' @export 
ejAttributes <- function(x){
	if(class(x) == "ejObject"){
		return(x$metadata)
	}
	if(class(x) == "ejRecord"){
		return(x$attributes)
	}
	if(class(x) == "ejEvent"){
		return(x$attributes)
	}
	stop("Unknown object type to extract attributes from.")
}

