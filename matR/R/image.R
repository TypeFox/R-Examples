
#-----------------------------------------------------------------------------------------
#  heatmap-dendrogram of biom objects.
#-----------------------------------------------------------------------------------------

image.biom <- function(
	x, ..., 
	map=NULL, 
	rows=TRUE,
	columns=TRUE, 
	rerender=NULL) {

	x <- x [rows, columns]
	arg <- list(...)

####  inherits (rerender, 'heatmap')
####  	components:		rowDendrogram, colDendrogram
####  	pass as:		Rowv, Colv
####  
####  inherits (rerender, 'dclust')
####  	compute:		as.dendrogram (rerender)		
####  	and pass as:	Rowv (or Colv...)
####  
####  inherits (rerender, 'dist')
####  	compute:		as.dendrogram (<hclustfun> (rerender))	(or, t(rerender))
####  	and pass as:	Rowv (or Colv...)

	Rowv <- Colv <- TRUE
	if (inherits (rerender, 'heatmap')) {
		Rowv <- rerender$rowDendrogram
		Colv <- rerender$colDendrogram
	} else if (inherits (rerender, 'dclust')) {
		stop("\'rerender\' is not yet implemented for class \'dclust\'")
	} else if (inherits (rerender, 'list')) {
		stop("\'rerender\' is not yet implemented for class \'list\'")
	} else if (inherits (rerender, 'dist')) {
		stop("\'rerender\' is not yet implemented for class \'dist\'")
	} else if (!is.null (rerender))
		stop("invalid \'rerender\' argument")

####  apply metadata references if present
####  apply metadata mapping

	arg$labRow <- subRow (arg$labRow, x)
	arg$labCol <- subColumn (arg$labCol, x)
	arg [names (map)] <- parMap (x, map, arg)

	par <- resolve (arg, list(
		x = as.matrix (x, TRUE),
		Rowv = Rowv,
		Colv = Colv,
		main = NULL,
		col = rgb (colorRamp(c("red", "green")) (seq(0, 1, length = 20)), maxColorValue = 255),

		colsep = 1:ncol(x),									# separate all columns
#		rowsep = 1:nrow(x),									# omit: this takes forever!
		sepwidth = 0.01,									# cell-separating gap
#		sepcolor =

		trace = "none",										# no midlines

#		mar, 
		margins = c(5,8),									# space for col/row labels

		labCol = colnames(x),								# column labels
		cexCol = 0.6,										# column label size (small)
		labRow = rownames(x),								# row labels
		cexRow = 0.6,										# row label size (small)

		key = FALSE,										# no key
#		xlab, ylab

#		lmat, lhei, lwid									# keep dendrograms small
		lhei = c(1,10),
		lwid = c(1,2)))
	yy <- suppressPackageStartupMessages (do.call (gplots::heatmap.2, par))

####  add class and 'call'

	yy$call <- match.call()
	class (yy) <- c("heatmap", "list")
	invisible (yy)
	}
