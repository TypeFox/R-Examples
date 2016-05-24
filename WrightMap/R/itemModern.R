itemModern <-
function(thr, yRange = NULL, axis.items = "Items", show.thr.sym = TRUE, thr.sym.cex = .8, thr.sym.lwd = 1, thr.sym.pch = 23, thr.sym.col.fg = rgb(0, 0, 0, 0.3), thr.sym.col.bg = rgb(0, 0, 0, 0.3), show.thr.lab = TRUE, thr.lab.pos = c(2, 4), thr.lab.text = NULL, thr.lab.col = "black", thr.lab.cex = .5, thr.lab.font = 2,label.items.rows = 1,label.items.srt = 0,label.items = NULL,label.items.cex = 0.6,label.items.ticks = TRUE,axis.logits = "Logits",show.axis.logits = "R", oma = c(0,0,0,3),cutpoints = NULL,...) {
	
	thr <- as.matrix(thr)

	nI <- dim(thr)[1]
	
	if(is.null(yRange))
		yRange <- c(min(thr, na.rm = TRUE),max(thr, na.rm = TRUE))


	if (is.null(thr.lab.text)) {
		if (!is.null(colnames(thr))) 
			thr.lab.text <- as.data.frame(matrix(rep(colnames(thr), each = nrow(thr)), nrow = nrow(thr)))
		else thr.lab.text = col(thr)
	}
	
	if (is.null(label.items)) {
		if (!is.null(rownames(thr))) 
			label.items <- rownames(thr)
		else label.items <- c(paste("Item", seq(1:nI)))
	}
	par(oma = oma)
	#par(mar = c(0,0,0,0))
	par(mgp = c(2, 0.2, 0))


	plot(seq(1:nI), rep(0, nI), type = "n", axes = FALSE, xlab = axis.items, ylab = "", ylim = yRange, xlim = c(0.5, nI + 0.5), 
		cex.lab = .8, font.lab = 3)

	
	if(!is.null(cutpoints)) {
		cutLines(cutpoints,...)
	}

	box(bty = "o")
	usr <- par("usr")
	par(mgp = c(3, 1, 0))
	
	if (show.axis.logits == "R" | show.axis.logits == TRUE) {
		axis(4, las = 1, cex.axis = 0.7, font.axis = 2)
		mtext(axis.logits, side = 4, line = 1.5, cex = 0.8, font = 3)
	} else if (show.axis.logits == "L") {
		axis(2, las = 1, cex.axis = 0.7, font.axis = 2)
		mtext(axis.logits, side = 2, line = 1.5, cex = 0.8, font = 3)
	}
	

#############
	if (show.thr.sym) {

		points(row(thr), thr, ylim = yRange, type = "p", cex = thr.sym.cex, lwd = thr.sym.lwd, pch = as.matrix(thr.sym.pch), 
			col = as.matrix(thr.sym.col.fg), bg = as.matrix(thr.sym.col.bg))
	}

	if (show.thr.lab) {
		if (show.thr.sym) {
			pos <- thr.lab.pos
			if(length(pos) != length(thr)) {
				pos <- matrix(rep(rep_len(pos, ncol(thr)), nI), byrow = TRUE, ncol = ncol(thr))
				pos <- t(sapply(1:nrow(thr), function(x) pos[x, rank(thr[x, ])]))
			}
			text(row(thr), thr, labels = as.matrix(thr.lab.text), col = as.matrix(thr.lab.col), pos = pos, cex = thr.lab.cex, 
				font = thr.lab.font)
		} else {
			text(row(thr), thr, labels = as.matrix(thr.lab.text), col = as.matrix(thr.lab.col), cex = thr.lab.cex, font = thr.lab.font)
		}
	}


##############

	#par(mgp = c(3, 1, 0))

	if (label.items.rows == 1) {
		if (label.items.srt != 0) {
			text.adj = c(1, 1)
		} else {
			text.adj = c(0.5, 2)
		}
		text(seq(1:nrow(thr)), y = par("usr")[3] - .2, labels = label.items, srt = label.items.srt, adj = text.adj, xpd = TRUE, cex = label.items.cex)

		if (label.items.ticks) {

			axis(1, at = 1:nI, labels = FALSE, line = NA, tcl = -0.35)

		}

	}

	if (label.items.rows == 2) {

		text(seq(from = 1, to = nrow(thr), by = 2), y = par("usr")[3], labels = label.items[seq(from = 1, to = nrow(thr), by = 2)], 
			adj = c(0.5, 2.4), xpd = TRUE, cex = label.items.cex)

		text(seq(from = 2, to = nrow(thr), by = 2), y = par("usr")[3], labels = label.items[seq(from = 2, to = nrow(thr), by = 2)], 
			adj = c(0.5, 4.0), xpd = TRUE, cex = label.items.cex)

		if (label.items.ticks == TRUE) {

			axis(1, at = seq(from = 1, to = nI, by = 2), labels = FALSE, line = NA, tcl = -0.35)
			axis(1, at = seq(from = 2, to = nI, by = 2), labels = FALSE, line = NA, tcl = -0.9)

		}

	}

	if (label.items.rows == 3) {

		text(seq(from = 1, to = nrow(thr), by = 3), y = par("usr")[3], labels = label.items[seq(from = 1, to = nrow(thr), by = 3)], 
			adj = c(0.5, 2.4), xpd = TRUE, cex = label.items.cex)


		text(seq(from = 2, to = nrow(thr), by = 3), y = par("usr")[3], labels = label.items[seq(from = 2, to = nrow(thr), by = 3)], 
			adj = c(0.5, 4.0), xpd = TRUE, cex = label.items.cex)


		text(seq(from = 3, to = nrow(thr), by = 3), y = par("usr")[3], labels = label.items[seq(from = 3, to = nrow(thr), by = 3)], 
			adj = c(0.5, 5.4), xpd = TRUE, cex = label.items.cex)

		if (label.items.ticks == TRUE) {

			axis(1, at = seq(from = 1, to = nI, by = 3), labels = FALSE, line = NA, tcl = -0.35)
			axis(1, at = seq(from = 2, to = nI, by = 3), labels = FALSE, line = NA, tcl = -0.9)
			axis(1, at = seq(from = 3, to = nI, by = 3), labels = FALSE, line = NA, tcl = -1.4)

		}

	}
	
	
}
