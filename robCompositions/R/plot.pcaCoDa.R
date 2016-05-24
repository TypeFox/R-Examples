#' Plot method
#' 
#' Provides robust compositional biplots.
#' 
#' The robust compositional biplot according to Aitchison and Greenacre (2002),
#' computed from resulting (robust) loadings and scores, is performed.
#' 
#' @param x object of class \sQuote{pcaCoDa}
#' @param y ...
#' @param \dots ...
#' @return The robust compositional biplot.
#' @author M. Templ, K. Hron
#' @seealso \code{\link{pcaCoDa}}
#' @references Aitchison, J. and Greenacre, M. (2002). Biplots of compositional
#' data. \emph{Applied Statistics}, \bold{51}, 375-392. \
#' 
#' Filzmoser, P., Hron, K., Reimann, C. (2009) Principal Component Analysis for
#' Compositional Data with Outliers. \emph{Environmetrics}, \bold{20} (6),
#' 621--632.
#' @keywords aplot
#' @export
#' @method plot pcaCoDa
#' @examples
#' 
#' data(expenditures)
#' p1 <- pcaCoDa(expenditures)
#' p1
#' plot(p1)
#' 
#' 
#' ## with labels for the scores:
#' data(arcticLake)
#' rownames(arcticLake) <- paste(sample(letters[1:26],nrow(arcticLake), replace=TRUE), 
#'                               1:nrow(arcticLake), sep="")
#' pc <- pcaCoDa(arcticLake, method="standard")
#' plot(pc, xlabs=rownames(arcticLake))
#' 
#' 
plot.pcaCoDa <- function(x, y, ...){
	## biplot
	#z <- list()
	#z$scores <- x$scores
	#z$loadings <- x$loadings
	beschx <- if(x$method == "robust") "PC1 (clr-robust)" else "PC1 (clr-standard)"
	beschy <- if(x$method == "robust") "PC2 (clr-robust)" else "PC2 (clr-standard)"
#	biplot(x$princompOutputClr, xlab=beschx, ylab=beschy)
	biplot(x$princompOutputClr, xlab=beschx, ylab=beschy, ...)
}
