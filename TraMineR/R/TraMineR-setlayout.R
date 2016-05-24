## =====================================
## Setting the layout for TraMineR plots
## =====================================

TraMineR.setlayout <- function(nplot, prows, pcols, withlegend, axes, 
	legend.prop=NA) {

	## Backward compatibility
	if (withlegend==TRUE) withlegend <- "auto"

	if (is.na(pcols)) pcols <- min(nplot,2)
	if (is.na(prows)) prows <- ceiling(nplot/pcols)

	## Defining initial layout matrix
	pheight <- 1
	pwidth <- 1
	widths <- rep(pwidth/pcols,pcols)
	heights <- rep(pheight/prows,prows)

	layrow <- prows
	laycol <- pcols
	laymat <- matrix(1:(layrow*laycol), nrow=layrow, ncol=laycol, byrow=TRUE)
	
	axisp <- 0

	legpos=NULL
	freecells <- (prows*pcols)-nplot

	## =========================
	## Positioning of the legend
	## =========================
	if (withlegend=="auto") {
		if (freecells==0) {
			if (is.na(legend.prop)) legend.prop <- 0.15
			layrow <- layrow+1

			pheight <- pheight-legend.prop
			heights <- rep(pheight/prows,prows)
			heights <- c(heights,legend.prop)

			widths <- rep(pwidth/laycol,laycol)

			legpos="bottom"

			## Adding one row in the layout matrix for the legend
			laymat <- rbind(laymat, rep(nplot+1,ncol(laymat)))
		}
		else {
			legpos="center"
			heights <- rep(pheight/prows,prows)
			widths <- rep(pwidth/laycol,laycol)
		}
	}
	else if (withlegend=="right") {
		if (is.na(legend.prop)) legend.prop <- 0.25
		laycol <- laycol+1
		pwidth <- pwidth-legend.prop
		legpos="center"
		widths <- rep(pwidth/pcols,pcols)
		widths <- c(widths, legend.prop)
		heights <- rep(pheight/prows,prows)

		## Adding one row in the layout matrix for the legend
		laymat <- cbind(laymat, rep(nplot+1,nrow(laymat)))
	}	

	## if (axes %in% c("all","bottom")) axisp <- 1

	## On which plots the axes will appear
	if (axes=="bottom") {
		for (nc in 1:ncol(laymat))
			axisp <- c(axisp, max(laymat[laymat[,nc]<=nplot,nc]))
	} 
	else if (axes=="all") axisp <- 1:nplot

	
	## Returning a list with layout settings
	laylist <- list(laymat=laymat, widths=widths, heights=heights, axisp=axisp, legpos=legpos)

	return(laylist)
} 
