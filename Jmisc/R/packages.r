#' load add-on packages. If the packages can not be found, install.packages is called.
#' @name packages
#' @aliases packages
#' @title load packages with auto-installation
#' @param x name of the packages
#' @param ... arguments to install.packages
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @seealso \link{require}  \link{install.packages}
#' @export
#' @examples     
#' \dontrun{
#' packages("foreach")
#' }
packages<-function(x,...){
	name_of_x <- 	gsub('"',"",deparse(substitute(x)))
	if (!require(name_of_x,character.only=TRUE)){
		install.packages(pkgs=name_of_x,...)
		require(name_of_x,character.only=TRUE)
	}
}