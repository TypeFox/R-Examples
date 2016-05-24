RECYCLEWARNING <- NULL
.onLoad <- function(libname, pkgname){
  RECYCLEWARNING <<- gettext(tryCatch( (1:2)+(1:3),warning=function(w) w$message ))
}

.onAttach <- function(libname, pkgname) {
  	pkgEnv = pos.to.env(match('package:SoundexBR', search()))
  	packageStartupMessage('')
}




#' @encoding latin1
#' @title ASCII Characters Table
#'
#' @description To detect ASCII characters, we may need to specify them literally. This function helps identifying what character is in ascii format and what is not.
#'
#' @param x A string whose characters is to be checked against.
#'
#' @return Returns \code{TRUE} if ASCII character and \code{FALSE} otherwise.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#'
#' @export
ascii.table <- function(x){
  charclass <- paste0("[^]"
, "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
, "0123456789"
, " !\"#$%&'()*+,./:;<=>?@[\\^_`{|}~-"
, "]")
  !grepl(charclass,x)
}