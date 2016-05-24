###################################################################################
#' Load a SFA object.
#'
#' @param filename		Load list \code{sfaList} from this file name
#'
#' @return list \code{sfaList}
#'
#' @references  \code{\link{sfaSave}} 
#' @keywords internal
#' @export
###################################################################################
sfaLoad<- function (filename){
	sfaList<-list()
	load(file=filename)
	return(sfaList)
}

###################################################################################
#' Save a SFA object.
#'
#' @param sfaList 		A list that contains all information about the handled sfa-structure
#' @param filename		Save list \code{sfaList} to this file
#'
#' @references  \code{\link{sfaLoad}} 
#' @keywords internal
#' @export
###################################################################################
sfaSave<- function (sfaList, filename){
	save(sfaList,file=filename)
}


###################################################################################
#' Load a GAUSS object.
#'
#' @param filename		Load list \code{gauss} from this file
#'
#' @return list \code{gauss}
#'
#' @references  \code{\link{gaussSave}} 
#' @keywords internal
#' @export
###################################################################################
gaussLoad<- function (filename){
	gauss<-list()
	load(file=filename)
	return(gauss)
}

###################################################################################
#' Save a GAUSS object.
#'
#' @param gauss 		  A list that contains all information about the handled gauss-structure
#' @param filename		Save list \code{gauss} to this file
#'
#' @references  \code{\link{gaussLoad}} 
#' @keywords internal
#' @export
###################################################################################
gaussSave<- function (gauss, filename){
	save(gauss,file=filename)
}