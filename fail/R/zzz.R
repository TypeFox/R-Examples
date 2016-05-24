#' @import checkmate
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom utils methods
#' @importFrom BBmisc stopf
#' @importFrom BBmisc warningf
#' @importFrom BBmisc collapse
#' @importFrom BBmisc %nin%
#' @importFrom BBmisc is.error
#' @importFrom BBmisc setClasses
#' @importFrom BBmisc names2
#' @importFrom BBmisc argsAsNamedList
#' @importFrom BBmisc vlapply

UNITCONVERT = c(1L, 1024L, 1048576L, 1073741824L)
names(UNITCONVERT) = c("b", "kB", "Mb", "Gb")
