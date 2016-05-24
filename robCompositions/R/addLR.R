#' @importFrom grDevices col2rgb
#' @importFrom grDevices dev.cur   
#' @importFrom grDevices dev.interactive
#' @importFrom grDevices gray
#' @importFrom grDevices grey
#' @importFrom grDevices rgb
#' @importFrom grDevices devAskNewPage
#' @importFrom graphics abline
#' @importFrom graphics identify
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics locator
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics points
#' @importFrom graphics segments
#' @importFrom graphics text
#' @importFrom stats aggregate  
#' @importFrom stats complete.cases  
#' @importFrom stats cov    
#' @importFrom stats cov.wt  
#' @importFrom stats density  
#' @importFrom stats dnorm  
#' @importFrom stats fitted  
#' @importFrom stats lm   
#' @importFrom stats lowess  
#' @importFrom stats mad 
#' @importFrom stats mahalanobis
#' @importFrom stats median
#' @importFrom stats model.matrix
#' @importFrom stats na.omit
#' @importFrom stats napredict
#' @importFrom stats pchisq
#' @importFrom stats pnorm
#' @importFrom stats princomp
#' @importFrom stats qchisq
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom stats residuals
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @importFrom stats terms
#' @importFrom stats var
#' @importFrom utils flush.console
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
NULL
 
#' Additive log-ratio transformation
#' 
#' The additive log-ratio transformation moves D-part compositional data from
#' the simplex into a (D-1)-dimensional real space.
#' 
#' The compositional parts are divided by the rationing part before the
#' logarithm is taken.
#' 
#' @param x D-part compositional data
#' @param ivar Rationing part
#' @return A list of class \dQuote{alr} which includes the following content:
#' \item{x.alr}{the transformed data} \item{varx}{the rationing variable}
#' \item{ivar}{the index of the rationing variable, indicating the column
#' number of the rationing variable in the data matrix \emph{x}}
#' \item{cnames}{the column names of \emph{x}} The additional information such
#' as \emph{cnames} or \emph{ivar} is usefull when a back-transformation is
#' applied on the \sQuote{same} data set.
#' @author Matthias Templ
#' @seealso \code{\link{addLRinv}}, \code{\link{isomLR}}
#' @references Aitchison, J. (1986) \emph{The Statistical Analysis of
#' Compositional Data} Monographs on Statistics and Applied Probability.
#' Chapman \& Hall Ltd., London (UK). 416p.
#' @keywords manip
#' @export
#' @examples
#' 
#' data(arcticLake)
#' x <- arcticLake
#' x.alr <- addLR(x, 2)
#' y <- addLRinv(x.alr)
#' ## This exactly fulfills:
#' addLRinv(addLR(x, 3))
#' data(expenditures)
#' x <- expenditures
#' y <- addLRinv(addLR(x, 5))
#' head(x)
#' head(y)
#' ## --> absolute values are preserved as well.
#' 
#' ## preserve only the ratios:
#' addLRinv(x.alr, ivar=2, useClassInfo=FALSE)
#' 
#' 
addLR <- function (x, ivar=ncol(x)){

	if(dim(x)[2] < 2) stop("data must be of dimension greater equal 2")
	x.alr <- log(x/x[, ivar])
	res <- list(x.alr=x.alr[,-ivar], 
			varx=x[,ivar], ivar=ivar, cnames=colnames(x))
	class(res) <- "alr"
	return(res)
}
