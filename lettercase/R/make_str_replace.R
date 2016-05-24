#' Make replace, delete and "is" functions for strings
#' 
#' Functions for building string functions for replacement, deletion and testing
#' 
#' @param pattern pattern to look for, as defined by a POSIX regular expression. 
#' See the "Extended Regular Expressions" section of regex for details. See 
#' fixed, ignore.case and perl for how to use other types of matching: fixed, 
#' case insensitive and perl-compatible expressions.
#' 
#' @param replacement string. References of the form \code{\\1}, \code{\\2} will be replaced 
#' with the contents of the respective matched group (created by ()) within the 
#' pattern.
#' 
#' These functions build functions that take a single string argument and return
#' a vector as a result.
#' 
#' \code{make_str_replace} builds a functions that replaces the strings 
#' according to the \code{pattern} and \code{replacement} arguments.
#' 
#' \code{make_str_delete} builds a functions that deletes the \code{pattern} 
#' from the string.
#' 
#' \code{make_str_is} builds a function that detects is the string is has a 
#' certain type of formatting.
#' 
#' @examples
#'   # -tk

   
make_str_replace <- function( pattern, replacement ) 
  function(string) {
    
    if (!is.atomic(string)) 
      stop("String must be an atomic vector", call. = FALSE)
    if (!is.character(string)) 
      string <- as.character(string)
  
    gsub( pattern, replacement, string, perl=TRUE) 
  }


#' @rdname make_str_replace
make_str_delete <- function( pattern )
  function(string) {
    
    if (!is.atomic(string)) 
      stop("String must be an atomic vector", call. = FALSE)
    if (!is.character(string)) 
      string <- as.character(string)
    
    gsub( pattern, '', string, perl=TRUE)
  }
