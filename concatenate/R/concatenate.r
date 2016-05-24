#' concatenate: Comma Concatenation
#'
#' Each function in concatenate returns a comma-separated string. (A length-one
#' character vector.) They can be used to construct human-friendly messages
#' whose elements aren't known in advance, like  calls to \code{message},
#' \code{warning} or \code{stop}, from clean code.
#' 
#' The workhorse function is \code{\link{cc}}. \code{\link{cn}} combines it with
#' grammatical number awareness, as in \code{\link{ngettext}}, and
#' \code{\link{sprintf}}-like substitution.
#'
#' @docType package
#' @name concatenate
NULL
