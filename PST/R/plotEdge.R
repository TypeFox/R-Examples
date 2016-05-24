## ==================================================
## Plotting the edges joining a node and its children
## ==================================================

plotEdge <- function(x1, x2, subtree, ePar, horiz = FALSE, gratio, max.level, cex, nc, cpal) {

	bx <- plotNodeLimit(x1, x2, subtree, max.level)
	xTop <- bx$x
	yTop <- subtree@order

	if (horiz) {
		X <- yTop
		Y <- xTop
		tmp <- yTop
		yTop <- xTop
		xTop <- tmp
	} else {
		Y <- yTop
		X <- xTop
	}

	stcol <- Xtract("stcol", ePar, default = cpal)
	type <- Xtract("type", ePar, default = "rectangle", 1)
	col <- Xtract("col", ePar, default = "grey", 1)
	pruned.col <- Xtract("pruned.col", ePar, default = col)
	lty <- Xtract("lty", ePar, default = par("lty"), 1)
	lwd <- Xtract("lwd", ePar, default = 3, 1)
	stem.height <- Xtract("stem.height", ePar, default=0.4)

	if (type == "rectangle") {
		## Bare verticale en dessous du rectangle
		if (horiz) {			
			segments(xTop, yTop, xTop+stem.height, yTop, col=col, lty=lty, lwd=lwd)
		} else {
			segments(X, yTop, X, yTop+stem.height, col=col, lty=lty, lwd=lwd)
		}
	}

	## Plotting edges to children
	children <- which.child(subtree)

	for (k in children) {
		child <- subtree[[k]]
		idx <- which(k==children)

		yBot <- child@order
		if (getOption("verbose")) { cat("    [-] plotEdge: to", child@path, "- L=", yBot, "\n") }
		xBot <- mean(bx$limit[idx:(idx + 1)])

		## REECRIRE UNE FONCTION is.leaf CORRRECTE PAR LA SUITE
		leaf <- length(which.child(child))==0

		i <- if (!leaf && child@order<max.level) 1 else 2

		## edge parameters
		c.col <- Xtract("c.col", ePar, default = "state", i)
		c.border <- Xtract("c.border", ePar, default = par("fg"), i)
		p.lwd <- Xtract("p.lwd", ePar, default = lwd, i)
		p.lty <- Xtract("p.lty", ePar, default = lty, i)
		## Color from circle to prob barplot
		ctp.col <- Xtract("ctp", ePar, default = "edge", i)

		if (type == "triangle") {
			if (horiz) {
				tmp <- xBot
				xBot <- yBot
				yBot <- tmp
				yNode.bottom <- yTop
			} 
			segments(xTop, yTop, xBot, yBot, col=col, lty=lty, lwd=lwd)
		} else {
			if (getOption("verbose")) { 
				cat("Child:", child@path, "yTop=", yTop, "xBot=", xBot, "xTop=", xTop)
			}

			## The color of the edge
			ecol <- if (ctp.col=="state" & k %in% names(stcol)) { stcol[k] } else { col }
			ccol <- if (c.col=="state" & k %in% names(stcol)) { stcol[k] } else { col }

			if (horiz) {
				## Bare horizontale du rateau
				segments(xTop+stem.height, xBot, xTop+stem.height, yTop, col=col, lty=lty, lwd=lwd)

				## Edge from the node to edge from parent
				segments(xTop+stem.height, xBot, xTop+1, xBot, col=ecol, lty=lty, lwd=lwd)
			} else {
				segments(xTop, yTop+stem.height, xBot, yTop+stem.height, col=col, lty=lty, lwd=lwd)

				## Edge from the node to edge from parent
				segments(xBot, yTop+stem.height, xBot, yTop+1, col=ecol, lty=lty, lwd=lwd)
			}
		}
	}
}
