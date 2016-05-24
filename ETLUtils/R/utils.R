###############################################################################
## Project: ETLUtils
## Content: Utilities like recoding, renaming
##
## Author: jan
## Creation date: Mar 25, 2012, 3:01:21 PM
## File: utils.R
###############################################################################


#' Recodes the values of a character vector
#'
#' Recodes the values of a character vector
#'
#' @param x character vector
#' @param from character vector with old values
#' @param to character vector with new values
#' @return 
#' x where from values are recoded to the supplied to values
#' @export
#' @seealso \code{\link{match}}
#' @examples
#' recoder(x=append(LETTERS, NA, 5), from = c("A","B"), to = c("a.123","b.123")) 
recoder <- function(x, from=c(), to=c()){
	missing.levels <- unique(x)
	missing.levels <- missing.levels[!missing.levels %in% from]
	if(length(missing.levels) > 0){
		from <- append(x=from, values=missing.levels)
		to <- append(x=to, values=missing.levels)
	}
	to[match(x, from)]
}




#' Performs NA replacement by last observation carried forward but adds 1 to the last observation carried forward. 
#' 
#' Performs NA replacement by last observation carried forward but adds 1 to the last observation carried forward. \cr
#'
#' @param x a numeric vector
#' @return a vector where NA's are replaced with the LOCF + 1 
#' @export
#' @seealso \code{\link[zoo]{na.locf}}
#' @examples
#' require(zoo)
#' x <- c(2,NA,NA,4,5,2,NA)
#' naLOCFPlusone(x)
naLOCFPlusone <- function(x){
	ix <- cumsum(is.na(x))
	zoo::na.locf(x) + ix - cummax(ix * !is.na(x))
}




#' Renames variables in a data frame.
#' 
#' Renames variables in a data frame. \cr
#'
#' @param x data frame to be modified.
#' @param from character vector containing the current names of each variable to be renamed.
#' @param to character vector containing the new names of each variable to be renamed.
#' @return The updated data frame x where the variables listed in from are 
#' renamed to the corresponding to column names. 
#' @export
#' @seealso \code{\link{colnames}, \link{recoder}}
#' @examples
#' x <- data.frame(x = 1:4, y = LETTERS[1:4])
#' renameColumns(x, from = c("x","y"), to = c("digits","letters"))
renameColumns <- function(x, from = "", to = ""){
	if(!"data.frame" %in% class(x)){
		stop("x should be of class data.frame")
	}
	if(length(from) != length(to)) {
		stop("from and to are not of equal length")
	}
	colnames(x) <- recoder(colnames(x), from = from, to = to)
	x
}



#' Put character vectors, columns of a data.frame or list elements as factor
#' 
#' Put character vectors, columns of a data.frame or list elements as factor if they are character strings
#' or optionally if they are logicals\cr
#'
#' @param x a character vector, a data.frame or a list
#' @param logicals logical indicating if logical vectors should also be converted to factors. Defaults to FALSE.
#' @param ... optional arguments passed on to the methods
#' @return The updated x vector/data.frame or list where the character vectors or optionally logical elements are 
#' converted to factors
#' @export
#' @seealso \code{\link{as.factor}}, \code{\link{factor}}
#' @examples
#' x <- data.frame(x = 1:4, y = LETTERS[1:4], b = c(TRUE, FALSE, NA, TRUE), stringsAsFactors=FALSE)
#' str(factorise(x))
#' str(factorise(x, logicals = TRUE))
#' str(factorise(list(a = LETTERS, b = 1:10, c = pi, d = list(x = x))))
factorise <- function(x, logicals=FALSE, ...){
  UseMethod(generic="factorise")
}
#' @rdname factorise
#' @export
factorise.default <- function(x, logicals=FALSE, ...){
  x
}
#' @rdname factorise
#' @export
factorise.character <- function(x, logicals=FALSE, ...){
  factor(x)
}
#' @rdname factorise
#' @export
factorise.data.frame <- function(x, logicals=FALSE, ...){
  if(logicals){
    tofactor <- sapply(x, FUN=function(x) inherits(x, "logical"))
    tofactor <- names(tofactor)[tofactor == TRUE]
    for(column in tofactor){
      x[[column]] <- factor(x[[column]])
    }
  }
  tofactor <- sapply(x, FUN=function(x) inherits(x, "character"))
  tofactor <- tofactor[tofactor == TRUE]
  if(length(tofactor) == 0){
    return(x)
  }
  x <- unclass(x)
  x <- as.data.frame(x, stringsAsFactors=TRUE)
  x
}
#' @rdname factorise
#' @export
factorise.list <- function(x, logicals=FALSE, ...){
  lapply(x, factorise, logicals = logicals)
}

