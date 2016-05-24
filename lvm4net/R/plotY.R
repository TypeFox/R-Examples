#' Plot the adjacency matrix of the network
#'
#' Function to plot the adjacency matrix of the network.
#'
#' @param Y list, or matrix containing a (\code{N} x \code{N}) binary adjacency matrix for each network view.
#' @param Ndata number of network views
#' @param EZ posterior mean latent positions
#' @param VZ posterior variance latent positions, if specified  draw ellipse
#' @param dimZ dimensions of Z to be plotted, default \code{dimZ = c(1, 2)}
#' @param labels text to be added in the plot representing the labels of each node. Default \code{labels = NULL}, no labels are shown 
#' @param colPl \code{col} for the points representing the nodes. Default \code{colPl = NULL}
#' @param colEll  \code{col} for the ellipses. Default \code{rgb(.6, .6 ,.6 , alpha=.1)}
#' @param LEVEL levels of confidence bounds shown when plotting the ellipses. Default \code{LEVEL = .95}
#' @param pchplot Default \code{pchplot = 20}
#' @param pchEll \code{pch} for the ellipses. Default \code{pchEll = 19}
#' @param pchPl \code{pch} for the points representing the nodes. Default \code{pchPl = 19}
#' @param cexPl \code{cex} for the points representing the nodes. Default \code{cexPl = 1.1}
#' @param arrowhead logical, if the arrowed are to be plotted. Default \code{arrowhead = FALSE}
#' @param curve curvature of edges. Default \code{curve = 0}
#' @param lwdLine lwd of edges. Default \code{lwdLine = .3}
#' @param xlim range for x 
#' @param ylim range for y
#' @param verbose if \code{verbose = TRUE} save the nodal positions 
#' @param ... Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}). 
#' @export
#' @examples
#' N <- 20
#' Y <- network(N, directed = FALSE)[,]
#' plotY(Y)
#' # Store the positions of nodes used to plot Y, in order to redraw the plot using 
#' # the same positions
#' z <- plotY(Y, verbose = TRUE)
#' plotY(Y, EZ = z)

plotY <- function(Y, Ndata = NULL, EZ = NULL, VZ = NULL, dimZ = c(1, 2), labels = NULL, colPl = 1, colEll = rgb(.6, .6 ,.6 , alpha=.1), LEVEL = .95, pchplot = 20, pchEll = 19, pchPl = 19, cexPl = 1.1, arrowhead = FALSE, curve = NULL,  lwdLine = .3,  xlim = NULL, ylim = NULL, verbose = FALSE, ...)
{	
	if(is.matrix(Y)) Y <- list(Y)
	
	stopifnot(is.list(Y), sapply(Y, is.adjacency))
	
	N <- nrow(Y[[1]])
	
	if(is.null(EZ)){ 
		EZ <- frEZ(Y[[1]], d = 2)
		
		stopifnot(dimZ %in% seq(1, 2) & length(dimZ == 2))
	} else {
		stopifnot(is.numeric(EZ) & is.matrix(EZ) & nrow(EZ) == N)
		stopifnot(dimZ %in% seq(1, ncol(EZ)) & length(dimZ == 2))
	}
		
	if(is.null(Ndata)){ 
		Ndata <- 1
		}
		
	if(length(Y) == 1 & Ndata > 1){
		for(i in 2:Ndata){
			Y[[i]] <- Y[[1]]
		}	
	}	
	
	stopifnot(length(Y) == Ndata)
		
		
	if(ncol(EZ) == 1) {

		EZ <- cbind(EZ, rep(0, N))
		
		if(is.null(curve)) curve <- .2
	
		}
		
	if(is.null(curve)) curve <- 0	
		
	Z <- EZ[,dimZ]
	
	if(is.null(xlim)) xlim <- range(Z[,1])
	if(is.null(ylim)) ylim <- range(Z[,2])

	plot(NA, xlim = xlim, ylim = ylim, xlab = expression(z[n1]), ylab = expression(z[n2]), 
       pch = pchplot, ...)
	
	if(!is.null(VZ)) 
	{
		for(n in 1:N)
		{
			coEl<-ellipse(VZ,centre = Z[n,],level=LEVEL,pch=pchEll)
			polygon(coEl[,1], coEl[,2], col = colEll,border=colEll)
		}
	}

	for(k in 1:Ndata){ 
		for(i in 1:(N-1)){ 
			for(j in 2:N){ 
				if(Y[[k]][i,j] == 1)
					network.arrow(Z[i,1], Z[i,2], Z[j,1], Z[j,2], col = rgb(0, 0, 0, alpha =.2), 
                        border = rgb(0, 0, 0, alpha =.2), lwd = lwdLine, 
                        arrowhead = arrowhead, curve = curve, ...)		
				}	
			}
		}
	
	if(!is.null(labels)){ 
		
		points(Z, pch = 3, cex= 1, col = colPl,...)
		text(Z, labels = labels, col = 1)
		
		} else {
			
				points(Z,pch = pchPl,col=colPl,...)
				points(Z,cex= cexPl,...)
		}
	
	if(verbose) return(Z)
	
}