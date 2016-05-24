## =======================================
## Plot method for objects of class pstree 
## =======================================

setMethod("plot", "PSTr", function (x, y=missing, max.level=NULL,
	nodePar = list(), edgePar = list(), 
	axis=FALSE, xlab = NA, ylab = if (axis) { "L" } else {NA}, horiz = FALSE,  
	xlim=NULL, ylim=NULL, 
	withlegend=TRUE, ltext=NULL, cex.legend=1, use.layout=withlegend!=FALSE, legend.prop=NA, ...) {

	cpal <- x@cpal

	## Margins
	Lmar <- if (axis) { 4 } else { 2 }
	if (horiz) { par(mar=c(Lmar,2,4,2)) } else { par(mar=c(2,Lmar,4,2)) }

	## Global setting for cex
	oolist <- list(...)
	cex <- if ("cex" %in% names(oolist)) { oolist[["cex"]] } else { 1 }

	## Axis are not plotted
	xaxt <- if ("xaxt" %in% names(oolist)) { oolist[["xaxt"]] } else { "n" }
	yaxt <- if ("yaxt" %in% names(oolist)) { oolist[["yaxt"]] } else { "n" }

	if (use.layout) {
		## Saving graphical parameters
		savepar <- par(no.readonly = TRUE)

		lout <- PST.setlayout(nplot=1, prows=NA, pcols=NA, withlegend, axes="all", legend.prop)
	  	layout(lout$laymat, heights=lout$heights, widths=lout$widths)

		legpos <- lout$legpos
	} else {
		legpos <- NULL
	}

	## List of segments
	seglist <- x@index

	## ORDER AT WHICH THE TREE STARTS
	k <- x@order
	stats <- summary(x, max.level=max.level, segmented=FALSE)
	
	if (missing(max.level) | is.null(max.level)) { max.level <- stats@depth }
	if (!"cpal" %in% names(nodePar)) { nodePar[["cpal"]] <- attr(x, "cpal") }

	hgt <- max.level
	mem.x <- stats@leaves
	pin <- par("pin")
	
	node.type <- Xtract("node.type", nodePar, default = "prob")
	node.size <- Xtract("node.size", nodePar, default = min(0.6, ((mem.x-1)/mem.x)-0.1))
	if (!"node.size" %in% names(nodePar)) { nodePar[["node.size"]] <- node.size }

	gratio <- Xtract("gratio", nodePar, default=min((((hgt-k)+1)/mem.x), 1))
	leave.lh <- Xtract("leave.lh", edgePar, default=0.1)
	leave.lw <- Xtract("leave.lw", edgePar, default=node.size)

	yTop <- k
	x1 <- 1
	x2 <- mem.x

	if (horiz) {
		xl. <- c(x1-((node.size/2)*gratio), x2 + ((node.size/2)*gratio))
		yl. <- c(k-(node.size/2), hgt+(node.size/2)+leave.lh)
	} else {
		ym <- if (node.type=="prob") { 0.5 } else { (node.size/2)*gratio }
		xl. <- c(x1 -(node.size/2), x2 + (node.size/2))
		yl. <- c(k-ym, hgt+((node.size/2)*gratio)+leave.lh)
	}
	yl. <- rev(yl.)

	## If horiz=TRUE, x and y are inverted
	if (horiz) {
		tmp <- xl.
		xl. <- yl.
		yl. <- rev(tmp)
		tmp <- xaxt
		xaxt <- yaxt
		yaxt <- tmp
		tmp <- xlab
		xlab <- ylab
		ylab <- tmp
	}
	if (missing(xlim) || is.null(xlim)) { xlim <- xl. }
	if (missing(ylim) || is.null(ylim)) { ylim <- yl. }

	plot(0, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, frame.plot=FALSE,
		ylab = ylab, xaxt=xaxt, yaxt=yaxt, ...)

	if (horiz) {
		## Size of one unit in inches (to correctly draw circles)  
		nc <- par("pin")[2]/par("usr")[3]
	} else {
		nc <- NULL
	}

	if (axis) {
		if (horiz) {
			axis(1, at=0:hgt)
		} else {
			axis(2, at=0:hgt)
		}		
	}

	plotTree(x1, x2, x, seglist, nPar = nodePar, ePar = edgePar, 
		horiz = horiz, gratio=gratio, max.level=max.level, cex=cex, nc=nc, cpal=cpal)

	## Plotting the legend
	if (!is.null(legpos)) {
		## Extracting some sequence characteristics
		if (is.null(ltext)) ltext <- x@labels

		PST.legend(legpos, ltext, cpal, cex=cex.legend)
	}

	## Restoring graphical parameters
	if (use.layout) {par(savepar)}
}
)

## 
Xtract <- function(nam, L, default, indx) {
	rep(if (nam %in% names(L)) L[[nam]] else default, length.out = indx)[indx]
}

##
asTxt <- function(x) {
	if (is.character(x) || is.expression(x) || is.null(x)) {x}
    		else {as.character(x)}
}





