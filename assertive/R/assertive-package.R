#' Readable check functions to ensure code integrity.
#'
#' \code{assertive} contains lots of \code{is_*} functions to check the
#' state of your variables, and \code{assert_*} functions to throw errors
#' if they aren't in the right form.
#'
#' When the package loads, it creates a global option 
#' \code{"assertive.severity"} that determines what happens when an 
#' \code{assert_*} function's input fails the condition.  By default, an error 
#' is thrown but it is possible to generate warnings or messages instead (see 
#' the examples).
#'
#' @docType package
#' @name assertive
#' @aliases assertive assertive-package
#' @examples
#' is_numeric(1:10)
#' assert_all_are_positive(1:10)
#' dont_stop(assert_is_scalar(runif(10)))
#' @author Richard Cotton \email{richierocks@@gmail.com}
NULL

#' Important changes to assertive
#' 
#' Changes since 0.3-0:
#' 
#' @section Virtuality:
#' The assertive package is now a virtual package: that is, it no longer 
#' contains its own functions, but instead rexports them from lower-level 
#' packages.
#' 
#' For interactive use, you can carry on using assertive as before.
#' 
#' For programmatic use, you can have more fine-grained control over what gets 
#' loaded by using the lower-level packages.
#' 
#' \code{assertive.base} contains the core functionality.
#' \code{assertive.properties} contains checks on properties of variables.
#' \code{assertive.types} contains checks on types of variables.
#' \code{assertive.numbers} contains checks for numbers.
#' \code{assertive.strings} contains checks for strings. 
#' \code{assertive.datetimes} contains checks for dates and times.
#' \code{assertive.files} contains checks for files and connections.
#' \code{assertive.sets} contains checks for sets.
#' \code{assertive.matrices} contains checks for matrices.
#' \code{assertive.models} contains checks for models.
#' \code{assertive.data} contains checks for complex data types.
#' \code{assertive.data.uk} contains checks for UK-specific complex data types.
#' \code{assertive.data.us} contains checks for US-specific complex data types. 
#' \code{assertive.reflection} contains checks on the state of R. 
#' \code{assertive.code} contains checks for code.
#' 
#' @section Translations:
#' The infrastructure for errors and warnings in multiple langauges is in place,
#' and translations are planned for future versions.  If you want to be a 
#' translator, email me at \email{richierocks@@gmail.com}.
#' @name changes
#' @docType package
NULL
