#' Centred log-ratio transformation
#' 
#' The cenLR transformation moves D-part compositional data from the simplex
#' into a D-dimensional real space.
#' 
#' Each composition is divided by the geometric mean of its parts before the
#' logarithm is taken.
#' 
#' @param x multivariate data ideally of class data.frame or matrix
#' @return The transformed data, including \item{x.clr}{clr transformed data}
#' \item{gm}{the geometric means of the original composition.}
#' @note The resulting transformed data set is singular by definition.
#' @author Matthias Templ
#' @seealso \code{\link{cenLRinv}}, \code{\link{addLR}}, \code{\link{isomLR}},
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
#' head(isomLRinv(eclr$x.clr))
#' 
cenLR <- function(x){
	#if(dim(x)[2] < 2) stop("data must be of dimension greater equal 2")
	if(dim(x)[2] == 1){
		res <- list(x.clr=x, gm=rep(1,dim(x)[1]))	    	
	} else{
		geometricmean <- function (x) {
			if (any(na.omit(x == 0)))
				0
			else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
		}
		gm <- apply(x, 1, geometricmean)
		x.clr <- log(x/gm)
		res <- list(x.clr=x.clr, 
				gm=gm
		)
	}
	class(res) <- "clr"
	return(res)  
}
