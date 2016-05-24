#' Inverse centred log-ratio transformation
#' 
#' Applies the inverse centred log-ratio transformation.
#' 
#' 
#' @param x an object of class \dQuote{clr}, \dQuote{data.frame} or
#' \dQuote{matrix}
#' @param useClassInfo if the object is of class \dQuote{clr}, the useClassInfo
#' is used to determine if the class information should be used. If yes, also
#' absolute values may be preserved.
#' @return the transformed data set.
#' @author Matthias Templ
#' @seealso \code{\link{cenLR}}, \code{\link{addLR}}, \code{\link{isomLR}},
#' \code{\link{addLRinv}}, \code{\link{isomLRinv}}
#' @references Aitchison, J. (1986) \emph{The Statistical Analysis of
#' Compositional Data} Monographs on Statistics and Applied Probability.
#' Chapman \& Hall Ltd., London (UK). 416p.
#' @keywords manip
#' @export
#' @examples
#' 
#' data(expenditures)
#' eclr <- cenLR(expenditures)
#' inveclr <- cenLRinv(eclr)
#' head(expenditures)
#' head(inveclr)
#' head(cenLRinv(eclr$x.clr))
#' 
cenLRinv <- function(x, useClassInfo=TRUE){
	
	if(class(x) != "clr" & useClassInfo == TRUE) warning("useClassInfo was set to FALSE, because x is not from class clr")
	if(!(class(x) %in% c("clr", "data.frame", "matrix"))) stop("class from x must be either clr, data.frame or matrix")
	xclr <- x$x.clr
	if(class(x) == "clr" & useClassInfo==TRUE){
		dat <- exp(xclr)  * x$gm
	} else if(class(x) != "clr" | useClassInfo==FALSE){
		dat <- exp(x)  		
	}
	return(dat)	
}
