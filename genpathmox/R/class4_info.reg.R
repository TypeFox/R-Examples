#' info.reg class
#'
#' info.pls is a S4 class that contains info on the variable and his levels that
#' provides the best binary split and the the the Fischers statitistcs: F-global,
#' F-coefficientes
#'
#' @name info.reg_class
#' @rdname info.reg_class

setClass("info.reg",representation(
variable		= "character",
modalidad		= "character",
fgstatistic   	= "numeric",
fpvalg         	= "numeric",
fcstatistic   	= "numeric",
fpvalc         	= "numeric",
candidates 		= "data.frame"
))
