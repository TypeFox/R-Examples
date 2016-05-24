#' Adjusting for original scale
#' 
#' Results from the model based iterative methods provides the results in
#' another scale (but the ratios are still the same). This function rescale the
#' output to the original scale.
#' 
#' It is self-explaining if you try the examples.
#' 
#' @param x object from class \sQuote{imp}
#' @return The object of class \sQuote{imp} but with the adjusted imputed data.
#' @author Matthias Templ
#' @export
#' @seealso \code{\link{impCoda}}
#' @references Hron, K. and Templ, M. and Filzmoser, P. (2010) Imputation of
#' missing values for compositional data using classical and robust methods
#' \emph{Computational Statistics and Data Analysis}, In Press, Corrected
#' Proof, ISSN: 0167-9473, DOI:10.1016/j.csda.2009.11.023
#' @keywords manip
#' @examples
#' 
#' data(expenditures)
#' x <- expenditures
#' x[1,3] <- x[2,4] <- x[3,3] <- x[3,4] <- NA
#' xi <- impCoda(x)
#' x
#' xi$xImp
#' adjust(xi)$xImp
#' 
adjust <- function(x){
	# x ... object from class "imp"
	if(class(x) != "imp") stop("object x must be from class imp")
	xneu=x$xImp
	s1 <- rowSums(x$xOrig, na.rm=TRUE)
	for(i in 1:nrow(x$xImp)){
		s <- sum(x$xImp[i, !x$wind[i,]])
		s2 <- sum(x$xImp[i, x$wind[i,]])
		fac <- s / (s + s2)
		s1[i] <-  s1[i] / fac
	}
	impS <- s1/rowSums(x$xImp)
	for(i in 1:ncol(x$xImp)){
		xneu[,i] <- x$xImp[,i] * impS
	}
	x$xImp <- xneu
	invisible(x)
}
