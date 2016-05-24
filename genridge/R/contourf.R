# Add facility to fill contours in contour() and to return the contourLines()

#  All arguments of contour.default commented out can be passed to contour() as ...
#  fill.col is a vector of fill colors corresponding to levels,
#    or a call to the color.palette() function.  
#    if not supplied, calculate a set of transparent colors 
#    based on col and fill.alpha
#  DONE:  add logic for fill.col=NULL to nothing extra, other than return contourLines
#  DONE:  make fill.alpha accept a more standard 0-1 numeric

contourf <- function(
		x = seq(0, 1, length.out = nrow(z)),
    y = seq(0, 1, length.out = ncol(z)),
    z,
    nlevels = 10, 
    levels = pretty(zlim, nlevels),
#    labels = NULL,
#    xlim = range(x, finite = TRUE),
#    ylim = range(y, finite = TRUE),
    zlim = range(z, finite = TRUE),
#    labcex = 0.6, drawlabels = TRUE, method = "flattest",
#    vfont, axes = TRUE, frame.plot = axes,
    col = par("fg"),
    color.palette = colorRampPalette(c("white", col)),
    fill.col = color.palette(nlevels+1),
    fill.alpha = 0.5,   # alpha transparency
#    lty = par("lty"), lwd = par("lwd"),
    add = FALSE, ...) {

	contour(x,y,z, nlevels=nlevels, levels=levels, zlim=zlim, col=col, add=add, ...)
	line.list <- contourLines(x, y, z, nlevels=nlevels, levels=levels)
	# contourLines returns a list of lists, each with components 
	# 'level', 'x', 'y'

	if (!is.null(fill.col)) {
		if(!is.na(fill.alpha)) {
			if (is.numeric(fill.alpha) && fill.alpha>=0 && fill.alpha<=1)
				fill.alpha <- as.hexmode(round(255 * fill.alpha))
			fill.col <- paste(fill.col, fill.alpha, sep="")
			}

		levs <- sapply(line.list, function(x) x[[1]])
		for (i in seq_along(line.list)) {
			clev <- which(levs[i] == unique(levs))
			polygon(line.list[[i]][2:3], col=fill.col[clev], border=NA)
		}
	}
	invisible(line.list)
}


