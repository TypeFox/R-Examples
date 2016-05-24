#' @title
#' Expect that the tested function is self-contained
#' 
#' @description
#'   Tests if a fuction is self-contained (i.e. do not use any global variables).
#' 
#' @param object 
#'   Function to test if it is self-contained.
#' @param label
#'   For full form, label of expected object used in error messages. 
#'   Useful to override default (deparsed expected expression) when doing 
#'   tests in a loop. For short cut form, object label. When NULL, computed from 
#'   deparsed object.
#' @param info 
#'   Extra information to be included in the message (useful when writing tests in loops).
#' 
#' @keywords internal
#' 
#' @export
expect_function_self_contained <- function(object, info = NULL, label = NULL) {
  lab <- make_label(object, label)

  global_vars <- codetools::findGlobals(object, merge = F)$variables

  if(length(global_vars)==0){
    testthat::succeed()
  } else {
    msg <- sprintf("%s contain global variable(s): %s.", lab, paste(global_vars, collapse = " "))
    testthat::fail(paste0(msg, info))
  }

  invisible(object)
}

#' @title Depricated function: expect_self_contained
#' 
#' @description Function has been depricated and will be removed. Please use \code{\link{expect_function_self_contained}} instead.
#' 
#' @keywords internal
#' @export
expect_self_contained <- function(object, info = NULL, label = NULL){
  .Deprecated("expect_function_self_contained")
  expect_function_self_contained(object, info, label)
}


#' @title
#' Expect that a given package is used
#' 
#' @description
#'   Tests that the following packages is used.
#' 
#' @param object
#'   Package to check for.
#' @param info 
#'   Extra information to be included in the message (useful when writing tests in loops).
#' 
#' @keywords internal
#' 
#' @export
expect_attached_package <- function(object, info = NULL){
  
  if(any(grepl(object, search()))){
    testthat::succeed()
  } else {
    msg <- sprintf("%s is not used.", object)
    testthat::fail(paste0(msg, info))
  }
  
  invisible(object)
}

#' @title Depricated function: expect_package
#' 
#' @description Function has been depricated and will be removed. Please use \code{\link{expect_attached_package}} instead.
#' 
#' @keywords internal
#' @export
expect_package <- function(object, info = NULL, label = NULL){
  .Deprecated("expect_attached_package")
  expect_attached_package(object, info)
}




#' @title
#' Expect function arguments
#' 
#' @description
#'  Test that an function object has a function with given arguments.
#' 
#' @param object
#'   Function to check the arguments of.
#' @param expected
#'   Expected arguments in function.
#' @param label
#'   For full form, label of expected object used in error messages. 
#'   Useful to override default (deparsed expected expression) when doing 
#'   tests in a loop. For short cut form, object label. When NULL, computed from 
#'   deparsed object.
#' @param info 
#'   Extra information to be included in the message (useful when writing tests in loops).
#' @param expected.label Equivalent of \code{label} for shortcut form.
#' 
#' @keywords internal
#' 
#' @export
expect_function_arguments <- function(object, expected, info = NULL, label = NULL, expected.label = NULL) {
  
  lab_obj <- make_label(object, label)
  lab_exp <- make_label(expected, expected.label)
  
  function_arguments <- names(formals(object))
  missing_arguments <- !function_arguments %in% expected
  extra_arguments <- !expected %in% function_arguments
  
  if(!(any(missing_arguments) | any(extra_arguments))){
    testthat::succeed()
  } else {
    msg <- sprintf("%s contain arguments: %s, not %s", 
                   lab_obj, 
                   paste(function_arguments, collapse = " "), 
                   lab_exp)
    testthat::fail(paste0(msg, info))
  }
  
  invisible(object)
}


#' @title
#' Expect function contain code
#' 
#' @description
#'  Test that a given code code exists in function
#' 
#' @param object
#'   Function to check for mandatory code
#' @param expected
#'   Expected arguments in function.
#' @param label
#'   For full form, label of expected object used in error messages. 
#'   Useful to override default (deparsed expected expression) when doing 
#'   tests in a loop. For short cut form, object label. When NULL, computed from 
#'   deparsed object.
#' @param info 
#'   Extra information to be included in the message (useful when writing tests in loops).
#' @param expected.label Equivalent of \code{label} for shortcut form.
#' 
#' @keywords internal
#' 
#' @export
expect_function_code <- 
  function(object, expected, info = NULL, label = NULL, expected.label = NULL) 
  {
    
    lab_obj <- make_label(object, label)
    lab_exp <- make_label(expected, expected.label)
    
    body <- as.character(body(object))
    
    if(any(grepl(x = body, pattern = expected))){
      testthat::succeed()
    } else {
      paste0("'", expected, "' not found in function body.")
      msg <- sprintf("%s not found in the body of %s", 
                     lab_exp, 
                     lab_obj)
      testthat::fail(paste0(msg, info))
    }
    
    invisible(object)
  }


# Functions taken from testthat package (that is not exported)

make_label <- function(object, label = NULL) {
  label %||% label(object)
}

label <- function(obj) {
  x <- lazyeval::lazy(obj)$expr
  
  if (is.character(x)) {
    encodeString(x, quote = '"')
  } else if (is.atomic(x)) {
    format(x)
  } else if (is.name(x)) {
    paste0("`", as.character(x), "`")
  } else {
    chr <- deparse(x)
    if (length(chr) > 1) {
      chr <- paste(deparse(as.call(list(x[[1]], quote(...)))), collapse = "\n")
    }
    chr
  }
}

`%||%` <- function(a, b) if (is.null(a)) b else a
