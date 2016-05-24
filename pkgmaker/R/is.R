# General test utility functions to check the type of objects
# 
# Author: Renaud Gaujoux
# Creation: 30 Apr 2012
###############################################################################
#' @include unitTests.R
NULL

#' Testing Object Type 
#' 
#' @name is_something
#' @rdname is_something
#' 
#' @return \code{TRUE} or \code{FALSE}
NULL

#' \code{is_NA} tests if a variable is exactly NA (logical, character, numeric or integer)
#' 
#' @param x an R object
#' @rdname is_something
#' @export
is_NA <- function(x){ 
	identical(x, NA) || identical(x, as.character(NA)) || identical(x, as.numeric(NA)) || identical(x, as.integer(NA))
}

#' \code{isFALSE} Tests if a variable is exactly FALSE.
#' 
#' @rdname is_something
#' @seealso \code{\link{isTRUE}}
#' @export
isFALSE <- function(x) identical(x, FALSE)

#' \code{isNumber} tests if a variable is a single number
#' 
#' @rdname is_something
#' @export
isNumber <- function(x){ 
	is.numeric(x) && length(x) == 1
}

#' \code{isReal} tests if a variable is a single real number
#' 
#' @rdname is_something
#' @export
isReal <- function(x){ 
	isNumber(x) && !is.integer(x)
}

#' \code{isInteger} tests if an object is a single integer
#' @rdname is_something
#' @export
isInteger <- function(x){ 
	is.integer(x) && length(x) == 1
}


#' \code{isString} tests if an object is a character string. 
#' 
#' @param y character string to compare with.
#' @param ignore.case logical that indicates if the comparison 
#' should be case sensistive.
#' 
#' @rdname is_something
#' @export
isString <- function(x, y, ignore.case=FALSE){
	if( res <- is.character(x) && length(x) == 1L ){
		if( !missing(y) ){
			if( !isString(y) ) stop("Invalid argument 'y': must be a string itself.")
			if( ignore.case ) {
				x <- toupper(x)
				y <- toupper(y)
			}
			res <-  x == y
		}
	}
	res
}

#' \code{is.dir} tests if a filename is a directory.
#' 
#' @rdname is_something
#' @export
is.dir <- function(x) file_test('-d', x)


#' \code{is.file} tests if a filename is a file.
#' 
#' @rdname is_something
#' @export
is.file <- function(x) file_test('-f', x)

#' \code{hasNames} tests if an object has names.
#' 
#' @param all logical that indicates if the object needs all names non empty
#' @rdname is_something
#' @export
hasNames <- function(x, all=FALSE){
	nm <- names(x)
	if( length(x) == 0L ) TRUE
	else !is.null(nm) && (!all || !is.element('', nm) )
}

unit.test(hasNames, {
			
			add_names <- function(x) setNames(x, letters[1:length(x)])
			# vector
			checkTrue(hasNames( add_names(1:10) ))
			checkTrue(hasNames( add_names(1:10) , all=TRUE))
			checkTrue(hasNames( c(add_names(1:10),11) ))
			checkTrue(!hasNames( c(add_names(1:10),11) , all=TRUE))
			# list
			checkTrue(hasNames( add_names(list(1,2,3)) ))
			checkTrue(hasNames( add_names(list(1,2,3)) , all=TRUE))
			checkTrue(hasNames( c(add_names(list(1,2,3)),11) ))
			checkTrue(!hasNames( c(add_names(list(1,2,3)),11) , all=TRUE))
			
		})
