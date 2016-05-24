#' Two dimensional plot of Latent Space Joint Model output
#'
#' Function to plot an object of class \code{'lsjm'}
#'
#' @param x object of class \code{'lsjm'}
#' @param Y list containing a (\code{N} x \code{N}) binary adjacency matrix for each network view.
#' @param drawCB logical if \code{drawCB = TRUE} draw confidence bounds 
#' @param dimZ dimensions of the latent variable to be plotted. Default \code{dimZ = c(1, 2)}
#' @param plotZtilde if TRUE do the plot for the last step of LSM
#' @param colPl \code{col} for the points representing the nodes. Default \code{colPl = NULL}
#' @param colEll  \code{col} for the ellipses. Default \code{rgb(.6, .6 ,.6 , alpha=.1)}
#' @param LEVEL levels of confidence bounds shown when plotting the ellipses. Default \code{LEVEL = .95}
#' @param pchplot Default \code{pchplot = 20}
#' @param pchEll \code{pch} for the ellipses. Default \code{pchEll = 19}
#' @param pchPl \code{pch} for the points representing the nodes. Default \code{pchPl = 19}
#' @param cexPl \code{cex} for the points representing the nodes. Default \code{cexPl = 1.1}
#' @param mainZtilde title for single network plots TRUE do the plot for the last step of LSM
#' @param arrowhead logical, if the arrowed are to be plotted. Default \code{arrowhead = FALSE}
#' @param curve curvature of edges. Default \code{curve = 0}
#' @param xlim range for x 
#' @param ylim range for y
#' @param main main title
#' @param ... Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}). 
#' @export
#' @examples
#'## Simulate Undirected Network
#'   N <- 20
#'   Ndata <- 2
#'    Y <- list()
#'    Y[[1]] <- network(N, directed = FALSE)[,]
#'    ### create a new view that is similar to the original
#'   for(nd in 2:Ndata){
#'     Y[[nd]] <- Y[[nd - 1]] - sample(c(-1, 0, 1), N * N, replace = TRUE, 
#'     prob = c(.05, .85, .1))
#'     Y[[nd]] <- 1 * (Y[[nd]]  > 0 )
#'   diag(Y[[nd]]) <- 0
#'    }
#'
#' par(mfrow = c(1, 2))
#' z <- plotY(Y[[1]], verbose = TRUE, main = 'Network 1')
#' plotY(Y[[2]], EZ = z, main = 'Network 2')
#' par(mfrow = c(1, 1))
#'
#' modLSJM <- lsjm(Y, D = 2) 
#' plot(modLSJM, Y, drawCB = TRUE)
#' plot(modLSJM, Y, drawCB = TRUE, plotZtilde = TRUE)

plot.lsjm <- function(x, Y, drawCB = FALSE, dimZ = c(1, 2), plotZtilde = FALSE, colPl = 1, 
                      colEll = rgb(.6, .6 ,.6 , alpha=.1), LEVEL = .95, pchplot = 20, 
                      pchEll = 19, pchPl = 19, cexPl = 1.1,  mainZtilde = NULL, 
                      arrowhead = FALSE, curve = NULL, xlim = NULL, ylim = NULL, 
                      main = NULL, ...)
{		
		stopifnot(inherits(x, 'lsjm'))
		stopifnot(is.logical(drawCB) & length(drawCB) == 1)
		stopifnot(dimZ %in% seq(1: ncol(x$EZ)) && length(dimZ == 2))
		stopifnot(is.logical(plotZtilde) & length(plotZtilde) == 1)
	
	if(plotZtilde){
		
		Ndata <- length(x$xiT)
		
		par(mfrow = c(1, Ndata))
		
		if(is.null(mainZtilde)){ 
			mainZtilde <- paste('Network', 1:Ndata)
			} else {
				if(length(mainZtilde == 1)){ 
					mainZtilde <- rep(mainZtilde, Ndata)
					} else {
						stopifnot(length(mainZtilde) == Ndata)
					}
			}
		
		for(nd in 1:Ndata){
		
			if(drawCB){ 
				VZ <- x$lsmVZ[[nd]]
			} else {
				VZ <- NULL
			}
	
		plotY(Y[[nd]], Ndata = 1, EZ = x$lsmEZ[[nd]], VZ = VZ, dimZ = dimZ, colPl = colPl, colEll = colEll, LEVEL = LEVEL, pchplot = pchplot, pchEll = pchEll, pchPl = pchPl, cexPl = cexPl, arrowhead = arrowhead, curve = curve,  xlim = xlim, ylim = ylim, main = mainZtilde[nd], ...)
		
		}
		
		par(mfrow = c(1, 1))
		
	} else {
	
		if(drawCB){ 
		VZ <- x$VZ
		} else {
		VZ <- NULL
		}
	
	plotY(Y, Ndata = length(x$xiT), EZ = x$EZ, VZ = VZ, dimZ = dimZ, colPl = colPl, colEll = colEll, LEVEL = LEVEL, pchplot = pchplot, pchEll = pchEll, pchPl = pchPl, cexPl = cexPl, arrowhead = arrowhead, curve = curve, xlim = xlim, ylim = ylim, main = main, ...)
	
	}
	
}