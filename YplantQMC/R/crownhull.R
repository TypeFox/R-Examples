#'Calculates and plots the convex hull around the plant crown
#'
#'@description This function finds the convex hull (and its surface area and volume)
#'spanning the leaves of the 3D plant, using all coordinates of the leaf edges.
#'The result is smallest set of x, y, z points that defines the convex hull,
#'that is, the polyhedral surface that contains all other points, and is
#'convex.
#'
#'The implementation is a wrapper for the 'convhulln' function in the from
#'package 'geometry'.
#'
#'Optionaly, uses the 'rgl' package (see \code{\link{plot3d}}), to add a plot
#'of the hull to the current (rgl) device.  Opens a new device if none is
#'currently open. Uses the non-visible function \code{triangles3d} to plot the
#'hull, see details there.
#'
#'The convex hull is calculated with the qhull algorithm, see
#'\code{\link{convhulln}} and references therein.
#'
#'@param xyz An object of class 'plant3d', or a matrix with three columns (xyz
#'coordinates).
#'@param plotit Logical. If FALSE, returns only volume and surface area.
#'@param alpha Transparency (0-1).
#'@return A list with components 'crownvolume' and 'crownsurface', giving the
#'volume and surface of the convex hull.
#'@author Remko Duursma
#'@references \url{www.qhull.org},
#'
#'Duursma, R.A., D.S. Falster, F. Valladares, F.J. Sterck, R.W. Pearcy, C.H.
#'Lusk, K.M. Sendall, M. Nordenstahl, N.C. Houter, B.J. Atwell, N. Kelly,
#'J.W.G. Kelly, M. Liberloo, D.T. Tissue, B.E. Medlyn and D.S. Ellsworth. 2012.
#'Light interception efficiency explained by two simple variables: a test using
#'a diversity of small- to medium-sized woody plants. New Phytologist.
#'193:397-408.
#'@keywords misc
#'@examples
#'
#'
#'# Toona example (plant included in package).
#'crownhull(toona)
#'
#'
#'# Some xyz data:
#'coords <- matrix(runif(300,0,1),ncol=3)
#'library(rgl)
#'plot3d(coords, col="blue", size=3, axes=FALSE, box=FALSE, xlab="", ylab="", zlab="")
#'crownhull(coords)
#'
#' @export
#' @importFrom rgl triangles3d
#' @importFrom geometry convhulln
crownhull <- function(xyz, plotit=TRUE, alpha=0.8){

	if(inherits(xyz, "plant3d"))
		xyz <- do.call("rbind", lapply(xyz$leaves, function(x)x$XYZ))

	# construct the hull (gives area and volume)
	if(is.list(xyz) && !is.data.frame(xyz))
		p <- as.matrix(do.call("rbind", xyz))
	else
		p <- as.matrix(xyz)
	
	ch <- convhulln(p, "FA")

	# construct the hull for plotting
	if(plotit){
		ch2 <- t(convhulln(p, "Qt"))
		triangles3d(p[ch2,1],p[ch2,2],p[ch2,3],col="forestgreen",alpha=alpha)
	}

return(list(crownvolume=ch$vol, crownsurface=ch$area))
}



