#' @noRd
#' @name precintcon.min
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.min 
#' @title Minimun 
#' @description Calculates the minimun precipitation and its related date into the serie. 
#' @usage precintcon.min(object) 
#' @param object a daily or monthly precipitation serie loaded vie the \code{\link{read.data}} function.
#' @return An array containing the following values:
#' \enumerate{
#' \item the minimun precipitation into the serie.
#' \item the date of the first ocorrence of the minimun precipitation into the serie.
#' }  
#' @seealso \code{\link{precintcon.limits.analysis}} 
#' @keywords precipitation limits min
precintcon.min <- function(object) {
	
	if (is.element("precintcon.daily", class(object))) {
		
		m <- min(as.matrix(object[,3:33]), na.rm=T)
		
		i <- which(object == m, arr.ind=T)
		
		d <- paste(i[1,2] - 2, object[i[1,1], 2], object[i[1,1], 1], sep="/")
		
		return(c(m, d))
		
	} else if (is.element("precintcon.monthly", class(object))) {
		
		m <- min(object[,3], na.rm=T)
		
		i <- which(object[,3] == m)
		
		d <- paste(object[i[1],2], object[i[1],1], sep="/")
		
		return(c(m, d))		
		
	} else
		stop("Invalid object type. It should be an precintcon.daily or precintcon.monthly object.")
	
};