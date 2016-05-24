#' @noRd
#' @name precintcon.classification
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.classification 
#' @title Classification of a precipitation serie
#' @description Classifies a daily precipitation serie based in an \code{interval}. 
#' @usage precintcon.classification(object, interval = 1) 
#' @param object is a daily precipitation serie. 
#' @param interval is an interval in millimeters for classifying.
#' @return A data.frame (precintcon.classified) containing the following variables:
#' \itemize{
#' \item \code{initial.class} is the initial value of the class.
#' \item \code{final.class} is the final value of the class.
#' \item \code{n} is the absolute freqnecy, i.e., the number of days in each class, 
#' i.e., days with precipitation amount between \code{initial.class} and \code{final.class}.
#' }
#' @seealso 
#' \code{\link{ci}}
#' \code{\link{read.data}}
#' \code{\link{as.daily}} 
#' \code{\link{as.fd}}
#' @keywords precipitation concentration index
precintcon.classification <- function(object, interval = 1) {
	
	if (is.element("precintcon.daily", class(object))) {
		
		v = as.vector(t(object[,3:33]))
		
		result <- data.frame()

		for (i in seq(0, max((v[!is.na(v)] + interval)), by=interval)) {
			
			initial <- if (i == 0) i + 0.1 else i
			
			final   <- (i + interval - 0.1)

			s <- sum(table(v[!is.na(v) & v >= initial & v <= final]))

			if (s > 0) result <- rbind(result, c(initial, final, s))
			   
		}
		
		if (nrow(result) > 0)
			colnames(result) <- c("initial.class", "final.class", "n")
		
		class(result) <- c("data.frame", "precintcon.classified")
		
		return(result)
		
	} else
		stop("precintcon.classification --> input data should be of class \"precintcon.daily\"")
}
