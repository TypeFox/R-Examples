#' Two dimensional plot of the Latent Space Model output
#'
#' Function to plot an object of class \code{'lsm'}
#'
#' @param x object of class \code{'lsm'}
#' @param Y (\code{N} x \code{N}) binary adjacency matrix
#' @param drawCB draw confidence bounds 
#' @param dimZ dimensions of the latent variable to be plotted. Default \code{dimZ = c(1, 2)}
#' @param colPl \code{col} for the points representing the nodes. Default \code{colPl = NULL}
#' @param colEll  \code{col} for the ellipses. Default \code{rgb(.6, .6 ,.6 , alpha=.1)}
#' @param LEVEL levels of confidence bounds shown when plotting the ellipses. Default \code{LEVEL = .95}
#' @param pchplot Default \code{pchplot = 20}
#' @param pchEll \code{pch} for the ellipses. Default \code{pchEll = 19}
#' @param pchPl \code{pch} for the points representing the nodes. Default \code{pchPl = 19}
#' @param cexPl \code{cex} for the points representing the nodes. Default \code{cexPl = 1.1}
#' @param arrowhead logical, if the arrowed are to be plotted. Default \code{arrowhead = FALSE}
#' @param curve curvature of edges. Default \code{curve = 0}
#' @param xlim range for x 
#' @param ylim range for y
#' @param ... Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}). 
#' @export
#' @examples
#' N <- 20
#' Y <- network(N, directed = FALSE)[,]
#'
#' modLSM <- lsm(Y, D = 2) 
#' plot(modLSM, Y)
#'
#' # Plot with 95% CB
#' plot(modLSM, Y, drawCB = TRUE)
#' # Plot with 99% CB
#' plot(modLSM, Y, drawCB = TRUE, LEVEL = .99)

plot.lsm <- function(x, Y, drawCB = FALSE, dimZ = c(1, 2), colPl = 1, 
                     colEll = rgb(.6, .6 ,.6 , alpha=.1), LEVEL = .95, 
                     pchplot = 20, pchEll = 19, pchPl = 19, cexPl = 1.1, 
                     arrowhead = FALSE, curve = NULL, xlim = NULL, ylim = NULL, ...)
{		
		stopifnot(inherits(x, 'lsm'))
		stopifnot(is.adjacency(Y))
		stopifnot(is.logical(drawCB) & length(drawCB) == 1)
		stopifnot(dimZ %in% seq(1: ncol(x$lsmEZ)) && length(dimZ == 2))
	
	if(drawCB){ 
		VZ <- x$lsmVZ
		} else {
		VZ <- NULL
		}
	
	plotY(Y, Ndata = 1, EZ = x$lsmEZ, VZ = VZ, dimZ = dimZ, colPl = colPl, colEll = colEll, LEVEL = LEVEL, pchplot = pchplot, pchEll = pchEll, pchPl = pchPl, cexPl = cexPl, arrowhead = arrowhead, curve = curve, xlim = xlim, ylim = ylim, ...)
	
}