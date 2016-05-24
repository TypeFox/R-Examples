# Mix image() and heatmap(), with multiple row annotation and distinct default values
# Author : Sylvain Mareschal <maressyl@gmail.com>
heat.map <- function(
	expr,
	side = NULL,
	cex.col = NA,
	cex.row = NA,
	mai.left = NA,
	mai.bottom = NA,
	mai.right = 0.1,
	mai.top = 0.1,
	side.height = 1,
	side.col = NULL,
	col.heatmap = heat(),
	zlim = "0 centered",   # or "range", or numeric(2)
	norm = c("rows", "columns", "none"),
	norm.robust = FALSE,
	customLayout = FALSE,
	getLayout = FALSE,
	font = c(1, 3)
	) {
	# Arg check
	norm <- match.arg(norm)
	
	# Layout parameters
	if(!isTRUE(customLayout)) {
		if(!is.null(side) && ncol(side) > 0) {
			mat <- matrix(c(1:2), ncol=1)
			heights <- c(lcm(ncol(side)*side.height), 1)
			widths <- 1L
		} else {
			mat <- 1L
			heights <- 1L
			widths <- 1L
		}
	} else {
		# No layout call
		mat <- as.integer(NA)
		heights <- as.integer(NA)
		widths <- as.integer(NA)
	}
	
	# Stop returning layout
	if(isTRUE(getLayout)) {
		return(list(mat=mat, heights=heights, widths=widths))
	}
	
	# Side color initialization
	if(is.null(side)) {
		# No side color
		side <- matrix(character(0), ncol=0, nrow=0)
		pal.side <- character(0)
	} else {
		# Check
		if(any(! rownames(expr) %in% rownames(side))) stop("All 'expr' row names must be in 'side' row names")
		
		# To reverted matrix
		side <- as.matrix(side)
		side <- side[ rownames(expr) , , drop=FALSE ]
		side <- side[ , ncol(side):1 , drop=FALSE ]
		
		# Values to color (ignore custom hexadecimal colors)
		val.side <- unique(as.character(side))
		val.side <- sort(val.side[ !is.na(val.side) ])
		val.side <- grep("^#([0-9A-Fa-f]{2}){3,4}$", val.side, invert=TRUE, value=TRUE)
		
		# Attribute colors to values
		if(length(val.side) > 0) {
			if(is.null(side.col)) {
				# Default palettes
				if(length(val.side) > 8) { pal.side <- rainbow(n=length(val.side), v=0.8)
				} else                   { pal.side <- c("#2a80b9", "#c1392b", "#f39c11", "#52be80", "#c185db", "#fee203", "#bec3c7")[1:length(val.side)]
				}
			} else {
				# Custom function
				pal.side <- side.col(length(val.side))
			}
			
			# Use value as name
			names(pal.side) <- val.side
		} else {
			# Empty legend (if only custom colors are provided)
			pal.side <- character(0)
		}
	}

	# Layout call
	if(!isTRUE(customLayout)) {
		layout(mat=mat, heights=heights, widths=widths)
		on.exit(layout(1))
	}
	
	# Centering and scaling output heatmap (heatmap() uses norm.robust=FALSE)
	if(isTRUE(norm.robust)) {
		center <- median
		scale <- mad
	} else {
		center <- mean
		scale <- sd
	}
	if(norm == "columns")     { expr <- (expr - apply(expr, 1, center, na.rm=TRUE)) / apply(expr, 1, scale, na.rm=TRUE)         # samples
	} else if(norm == "rows") { expr <- t((t(expr) - apply(expr, 2, center, na.rm=TRUE)) / apply(expr, 2, scale, na.rm=TRUE))   # genes
	}
	
	# Symmetrical palette around 0
	if(identical(zlim, "0 centered")) {
		zlim <- max(abs(min(expr, na.rm=TRUE)), max(expr, na.rm=TRUE), na.rm=TRUE)
		zlim <- c(-zlim, zlim)
	} else if(identical(zlim, "range")) {
		zlim <- range(expr, na.rm=TRUE)
	}
	
	# Apply ceiling
	expr[ expr < zlim[1] ] <- zlim[1]
	expr[ expr > zlim[2] ] <- zlim[2]
	
	# Evolutive cex (from heatmap())
	if(is.na(cex.col)) cex.col <- 0.2 + 1 / log10(nrow(expr))
	if(is.na(cex.row)) cex.row <- 0.2 + 1 / log10(ncol(expr))
	
	# Evolutive margins
	if(is.na(mai.left))   mai.left   <- max(strwidth(colnames(expr), units="inches", cex=cex.row)) + par("cin")[2]
	if(is.na(mai.bottom)) mai.bottom <- max(strwidth(rownames(expr), units="inches", cex=cex.col)) + par("cin")[2]
	
	# Side color plot (top)
	if(ncol(side) > 0) {
		par(mai=c(0, mai.left, mai.top, mai.right))
		plot(x=NA, y=NA, xlim=c(0.5, nrow(expr)+0.5), ylim=c(0, ncol(side)), xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab="")
		for(k in 1:ncol(side)) {
			# Annotation colors
			col <- side[,k]
			isCustom <- grepl("^#([0-9A-Fa-f]{2}){3,4}$", col)
			col[!isCustom] <- pal.side[ col[!isCustom] ]
			
			# Draws annotation boxes
			rect(xleft=(1:nrow(expr))-0.5, xright=(1:nrow(expr))+0.5, ybottom=k-1L, ytop=k, col=col)
			
			# Add annotation title
			if(!is.null(colnames(side))) mtext(side=2, at=k-0.5, text=colnames(side)[k], las=2, line=1)
		}
	}
	
	# Heatmap (bottom)
	par(mai=c(mai.bottom, mai.left, mai.top, mai.right))
	image(expr, xaxt="n", yaxt="n", col=col.heatmap, zlim=zlim)
	axis(side=1, at=(1:nrow(expr) - 1L) / (nrow(expr) - 1L), labels=rownames(expr), las=2, cex.axis=cex.col, tick=FALSE, line=-0.5, font=font[1])
	axis(side=2, at=(1:ncol(expr) - 1L) / (ncol(expr) - 1L), labels=colnames(expr), las=2, cex.axis=cex.row, tick=FALSE, line=-0.5, font=font[2])
	box()
	
	# Invisibly return parameters for heatScale()
	invisible(
		list(
			zlim = zlim,
			col.heatmap = col.heatmap,
			legend = pal.side,
			cex.col = cex.col,
			cex.row = cex.row,
			mai.left = mai.left,
			mai.bottom = mai.bottom
		)
	)
}

