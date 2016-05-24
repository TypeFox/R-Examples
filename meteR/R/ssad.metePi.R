#' @title Species Spatial Abundance Distribution
#'
# @description 
#'
# @details
#' 
#' 
#' @param x An objects of class meteSSF; i.e. the spatial structure function \eqn{\Pi(n)}
#' 
#' @export
#' 
#' @examples
#' data(anbo)
#' pi1 <- meteSSF(anbo$spp, 'crcr', anbo$count, row=anbo$row, col=anbo$col, A=1, A0=16)
#' plot(ssad(pi1))

# @return list
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @seealso sad.mete, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

ssad <- function(x) {
  UseMethod('ssad')
}

#' @rdname ssad
# @method ssad meteSSF
# @S3method ssad meteSSF
#' @export 

ssad.meteSSF <- function(x) {
	dat <- x$data$n
	
	if(is.null(dat)) {
		X <- NULL
	} else {
		X <- sort(dat, decreasing=TRUE)
	}
	
	this.eq <- function(n, log=FALSE) {
		out <- metePi(n, x$La, x$state.var['n0'])
		if(log) out <- log(out)
		
		return(out)
	}
	
	FUN <- distr::DiscreteDistribution(supp=0:x$state.var['n0'],
	                                   prob=this.eq(0:x$state.var['n0']))
	
	out <- list(type='ssad', data=X,
	            d=this.eq, p=FUN@p, q=FUN@q, r=FUN@r,
	            state.var=x$state.var, La=x$La)
	
	class(out) <- c('ssad', 'meteDist')
	
	return(out)
}
