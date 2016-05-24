# Logging system
# 
# Author: Renaud Gaujoux
# Creation: 25 Apr 2012
###############################################################################

#' @include utils.R
NULL

#' Internal verbosity option
#' @param val logical that sets the verbosity level.
#' @return the old verbose level   
#' @keywords internal
lverbose <- local({
			.val <- NA
			function(val){
				if( missing(val) ) return(.val)
				oval <- .val
				.val <<- val
				invisible(oval)
			}
		})

#' Tells if all verbose messages should be 
#' @rdname lverbose
lsilent <- function(){
	l <- lverbose()
	is.na(l) || l == 0L
}
#' Tells if verbosity is on.
#' @rdname lverbose
is.verbose <- function(){
	l <- lverbose()
	!is.na(l) && l >= 0L
}

#' Prints out a message (on sdtout) if verbose mode is on.
#' 
#' @param ... arguments passed to \code{...} \code{\link{cat}}
#' @param appendLF logical indicating if an endline character should be appended 
#' at the end of the message. Passed to \code{\link{cat}}.
#' @rdname lverbose
#'  
vmessage <- function(...){
	lmessage(..., level=1L)
}
#' Prints out a message (on sdtout) if the verbosity level is greater than a 
#' given value. 
#' 
#' @param level verbosity level threshold (numeric value) above which the 
#' message should be printed out. 
#' This threshold is compared with the current verbosity level as returned by 
#' \code{lverbose}.
#' @param ... arguments passed to \code{...} \code{\link{lmessage}} or \code{\link{cat}} 
#' @param appendLF logical indicating if an endline character should be appended 
#' at the end of the message. Passed to \code{\link{cat}}.
#' 
#' @rdname lverbose
#' 
lmessage <- function(..., level=1L, appendLF=TRUE){
	l <- lverbose()
	if( !is.na(l) && l >= level ) cat(..., if(appendLF) "\n", sep='')
}


