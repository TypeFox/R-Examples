#'Plots a Yplant leaf file (a file with extension \code{.l} or \code{.lf}) in
#'2D.
#'
#'Produces a plot of a Yplant leaf file, read in using \code{\link{readl}}.
#'
#'
#'@aliases plot.leaffile plot3d.leaffile
#'@param x Object of class 'leaffile'.
#'@param nleaf Which leaf to plot in the leaf file (if more than one leaf
#'available in the file).
#'@param edgepoints Logical. If TRUE, plots dots on the leaf edge coordinates.
#'@param edgecex If edgepoint=TRUE, cex (i.e. size) of the leaf edge dots.
#'@param \dots Further parameters passed to \code{plot.default}.
#'@author Remko Duursma
#'@seealso \code{\link{readl}}
#'@keywords misc
#'@examples
#'
#'
#'\dontrun{
#'
#'# Read and plot a leaf in one go, select a leaf from a menu.
#'plot(readl())
#'
#'# Make a pdf of all leaf files in the current working directory:
#'leaffiles <- list.files(pattern="\\.l$", ignore.case=TRUE)
#'pdf("Leaf files.pdf", onefile=TRUE)
#'for(i in 1:length(leaffiles))plot(readl(leaffiles[i]))
#'dev.off()
#'
#'
#'}
#'
#'@export
#'@method plot leaffile
#'@S3method plot leaffile
plot.leaffile <- function(x, nleaf=1, edgepoints=TRUE, edgecex=0.8, ...){
	xyz <- x[[nleaf]]$XYZ
	leaflen <- max(xyz[,"Y"]) - min(xyz[,"Y"])
	leafwid <- max(xyz[,"X"]) - min(xyz[,"X"])
	span <- max(leaflen, leafwid)
	neglen <- min(xyz[,"Y"])
	Ylim <- c(neglen, span-neglen)
	Xlim <- c(-span/2, span/2)
	
	par(pty='s')
	plot(xyz[,"X"], xyz[,"Y"], type='n', pch=19, xlim=Xlim, ylim=Ylim,
	xlab="Width (mm)", ylab="Length (mm)",...)
	polygon(xyz[,"X"], xyz[,"Y"], col="forestgreen")
	if(edgepoints)points(xyz[,"X"], xyz[,"Y"], pch=19, cex=edgecex)
}



