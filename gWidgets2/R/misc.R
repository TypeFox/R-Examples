## miscellaneous functions

##' @import methods
##' @import digest
NULL

## merge two lists
## 
## @param x a list
## @param y a list
## @param ... passes along
## @note not exported, appears in some other packages. Call via \code{:::}.
merge.list <- function(x, y, ...) {
  args <- list(...)
  overwrite <- getWithDefault(args$overwrite, TRUE)
  
  if(missing(y) || is.null(y))
    return(x)
  nms <- names(y)
  for(i in seq_along(y)) {
    nm <- nms[i]                        # possibly NULL
    if(is.null(nm)) {
      x[[length(x) + 1]] <- y[[i]]
    } else if(overwrite || !(nm %in% names(x))) {
      x[[nm]] <- y[[i]]
    }
  }
  x
}

##' Flatten a nested list
##'
##' @author Tommy (http://stackoverflow.com/questions/8139677/how-to-flatten-a-list-to-a-list-without-coercion)
##' @param x a list
flatten <- function(x) {
  repeat {
    if(!any(vapply(x, is.list, logical(1)))) return(x)
    x <- Reduce(c, x)
  }
}


##' get index of element of list
##'
##' Like match, but works with a list
##' @param lst a list to search through
##' @param ele element of list
##' @return returns index of element or integer(0)
get_index_in_list <- function(lst, ele) {
  n <- seq_along(lst)
  n[sapply(lst, function(i) identical(i, ele))]
}

##' Return x unless NULL, NA, length 0, ..., in which case we give default
##'
##' @param x value to return or its default
##' @param default default value
##' @return x or default
##' @export
getWithDefault <- function(x, default) {
  if(is_empty(x))
    default
  else
    x
}

##' is value missing, null, 0-length or NA length 1
##'
##' @param x object to test
##' @return logical
##' @export
is_empty <- function(x) {
  missing(x) ||
  is.null(x) ||
  length(x) == 0 ||
  (is.atomic(x) && length(x) == 1 && is.na(x))
}

##' Functions to message something needs doing. Easy to search for
##'
##' @param msg optional message to emit
##' @export
XXX <- function(msg) {
  if(!missing(msg))
    message(msg)
}

##' Method to send message if any depreacted arguments are being used
##'
##' Many arguments were deprecated due to various reasons. This is meant to ease porting of code.
##' @param deprecated_args named list of deprecated args
##' @param ... named avlues
check_deprecated <- function(deprecated_args=list(), ...) {
  if(!length(deprecated_args))
    return()

  args <- list(...)
  sapply(names(args), function(i) {
    if(!is.null(tmp <- deprecated_args[[i, exact=TRUE]]))
      message(sprintf("Argument %s has been deprecated:\n\t o ",
                      i, paste(tmp, collapse="\n\t")))
  })
}

##' check that toolkit object return the right class
##' 
##' The S3 dispatch assumes naming conventions in the class names. This offers some check.
##' @param obj object with expected return class
##' @param ret_class character string of class expected
##' @return throws error if a mismatch
check_return_class <- function(obj, ret_class) {
  if(!any(sapply(ret_class, is, object=obj)))
    stop(sprintf("Expecting toolkit object of class (or subclass) %s. Got one of class %s",
                 paste(ret_class, collapse="; "), class(obj)[1]))
}
  
##' blurb about installation
##'
##' put in so can be updated easily
installing_gWidgets_toolkits <- function() {
  file <- system.file("install/installing_toolkits.txt", package="gWidgets2")
  tmp <- readLines(file)
  cat(paste(tmp, "\n"))
}


##' Return logical indicating if we are on a macintosh machine
##'
##' @return logical
##' @export
is_MacOSX <- function() {
  grepl("darwin", R.Version()$os)
}

##' Return logical indicating if we are on a Windows machine
##'
##' @return logical
##' @export
is_Windows <- function() {
  grepl("Windows", R.Version()$os)
}
  
## some special class unions so we can have easier to deal with default
setClassUnion("IntegerOrNULL", c("integer", "NULL"))
setClassUnion("CharacterOrNULL", c("character", "NULL"))
setClassUnion("LogicalOrNULL", c("logical", "NULL"))
setClassUnion("LogicalCharacterOrNULL", c("logical", "character", "NULL"))
setClassUnion("FunctionOrNULL", c("function", "NULL"))
              
